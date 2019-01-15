function [analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams] = fit_IAMP(analysisParams, fullCleanData)
% Takes in the clean time series data and the analysis params and fits the IAMP model.
%
% Syntax:
%   [analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams] = fit_IAMP(analysisParams, fullCleanData);
%
% Description:
%    This function takes in the clean time series data and the analysis params
%    and fits the IMAP model. This function builds a stimulus design matirx
%    based on the analysisParams (from each run of the experiemnt) and run the
%    IAMP model on the cleaned and trial sorted data. 
%
% Inputs:
%    analysisParams             - Struct of important information for the
%                                 analysis 
%    fullCleanData              - The cleaned time course
%
% Outputs:
%    analysisParams             - Returns analysisParams with any updates
%    iampTimeCoursePacketPocket - Cell array of IAMP packets for each run 
%    iampOBJ                    - The IAMP object
%    iampParams                 - Cell array of IAMP parameter fits for each run 
%
% Optional key/value pairs:
%    none

% MAB 09/09/18
% MAB 01/06/19 -- changed from runIAMP_QCM to fit_IAMP and removed QCM


analysisParams.numSessions = length(analysisParams.sessionFolderName);

for sessionNum = 1:analysisParams.numSessions
    
    % Gets the path to a text file that contains the mat file names needed
    % to get the trail order information for each run.
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'melaDataPath'), analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderFiles = textFile2cell(trialOrderFile);
    
    % Get the Directions of each session. This requires analysisParams.directionCoding to be organized
    % such that the directions are grouped per session and these groups are in the same order as the
    % sessions order
    
    % get then number of direction over all sessiosn
    analysisParams.numDirections = size(unique(analysisParams.directionCoding','rows')',2);
    
    % Get the directions and contrasts for corresponding session
    if analysisParams.numDirPerSession == analysisParams.numDirections
        directionCoding  = analysisParams.directionCoding;
        maxContrast = analysisParams.maxContrastPerDir;
    elseif analysisParams.numDirPerSession < size(unique(analysisParams.directionCoding','rows')',2)
        sPos = 1+ analysisParams.numDirPerSession*(sessionNum-1);
        ePos = (1+ analysisParams.numDirPerSession*(sessionNum-1)) + (analysisParams.numDirPerSession-1);
        directionCoding = analysisParams.directionCoding(:,sPos:ePos);
        maxContrast = analysisParams.maxContrastPerDir(sPos:ePos);
    else
        error('number of directions per session is greater than the number of total direction')
    end
    
    %% Construct the model object
    iampOBJ = tfeIAMP('verbosity','none');
    
    %% Create a cell of stimulusStruct (one struct per run)
    for jj = 1:analysisParams.numAcquisitions
        
        % identify the data param file
        dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
        
        % We are about to load the data param file. First silence the warning
        % for EnumerableClassNotFound. Save the warning state.
        warningState = warning();
        warning('off','MATLAB:class:EnumerableClassNotFound')
        
        % Load and process the data param file
        load(dataParamFile);
        expParams = getExpParams(dataParamFile,analysisParams.TR,'hrfOffset', false, 'stripInitialTRs', false);
        
        % restore warning state
        warning(warningState);
        
        % make timebase
        totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
        deltaT = analysisParams.TR*1000;
        stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
        responseStruct.timebase = stimulusStruct.timebase;
        
        % make stimulus values for IAMP
        % Stim coding: 80% = 1, 40% = 2, 20% = 3, 10% = 4, 5% = 5, 0% = 6;
        stimulusStruct.values =  createRegressors(expParams,analysisParams.baselineCondNum,totalTime,deltaT);
        
        % make stimulus values for QCM
        contrastCoding = [analysisParams.contrastCoding, 0];
        LMSContrastMat = LMSContrastValuesFromParams(expParams,contrastCoding,directionCoding,maxContrast,totalTime,deltaT);
        directionPrecision = 4;
        indDirectionDirections = round(directionCoding(1:analysisParams.theDimension,:),directionPrecision);
        LMSContrastMat(3,:) = [];
        [stimDirections,stimContrasts] = tfeQCMStimuliToDirectionsContrasts(LMSContrastMat, ...
            'zeroContrastDirection',indDirectionDirections(:,1),'precision',directionPrecision);
        
        % Set the number of instances.
        clear defaultParamsInfo
        defaultParamsInfo.nInstances = size(stimulusStruct.values,1);
        
        % Get the kernel
        kernelStruct = generateHRFKernel(6,12,10,stimulusStruct.timebase);
        
        % Take the median across voxels
        responses = median(fullCleanData(:,:,(jj+((sessionNum-1)*10))),1)';
        
        %%  Make the IAMP packet
        % the response
        thePacket.response.values   = responses';
        thePacket.response.timebase = stimulusStruct.timebase;
        % the stimulus
        thePacket.stimulus.timebase = stimulusStruct.timebase;
        thePacket.stimulus.values   = stimulusStruct.values;
        % the kernel
        thePacket.kernel = kernelStruct;
        % the meta data (this is the constrast and directions)
        thePacket.metaData.stimDirections = stimDirections;
        thePacket.metaData.stimContrasts  = stimContrasts;
        thePacket.metaData.lmsContrast    = LMSContrastMat;
        
        % Perform the fit
        [paramsFit,fVal,IAMPResponses] = ...
            iampOBJ.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, ...
            'searchMethod','linearRegression');
        
        iampParams{sessionNum,jj} = paramsFit;
        iampTimeCoursePacketPocket{sessionNum,jj} = thePacket;

    end
    
end

end


