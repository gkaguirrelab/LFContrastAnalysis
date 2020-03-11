function [analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData, varargin)
% Takes in the clean time series data and the analysis params and fits the IAMP model.
%
% Syntax:
%    analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams,
%                                                  fullCleanData,varargin);
%
% Description:
%    This function takes in the clean time series data and the analysis params
%    and creates a packet per run. The output is a cell array of packets
%    the length of the number of runs. 
%
% Inputs:
%    analysisParams             - Struct of important information for the
%                                 analysis
%    fullCleanData              - The cleaned time course
%
% Outputs:
%    analysisParams             - Returns analysisParams with any updates
%    iampTimeCoursePacketPocket - Cell array of packets for each run
%   
% Optional key/value pairs:
%    highpass                   - use a high pass filter on the data

% MAB 03/10/20 Wrote it.

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addRequired('fullCleanData',@isnumeric);
p.addParameter('plotColor',[],@isvector);
p.addParameter('highpass',false,@islogical)

p.parse(analysisParams,fullCleanData,varargin{:});

% Get the number of sessions
analysisParams.numSessions = length(analysisParams.sessionFolderName);

% Loop over sessions
for sessionNum = 1:analysisParams.numSessions
    
    % Gets the path to a text file that contains the mat file names needed
    % to get the trail order information for each run.
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
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
        regMat=  createRegressors(expParams,analysisParams.baselineCondNum,totalTime,deltaT);
        
        % this converts the 21xt regression matrix to 41xt martix to
        % code for the 40 total conditions + baseline
        if sessionNum == 1
            stimulusStruct.values = [regMat(1:end-1,:);zeros(size(regMat(1:end-1,:)));regMat(end,:)];
        elseif sessionNum == 2
            stimulusStruct.values = [zeros(size(regMat(1:end-1,:)));regMat];
        else
            error('A thrid session is not coded for')
        end
        
        % make stimulus values for QCM
        contrastCoding = [analysisParams.contrastCoding, 0];
        LMSContrastMat = LMSContrastValuesFromParams(expParams,contrastCoding,directionCoding,maxContrast,totalTime,deltaT);
        directionPrecision = 4;
        indDirectionDirections = round(directionCoding(1:analysisParams.theDimension,:),directionPrecision);
        LMSContrastMat(3,:) = [];
        [stimDirections,stimContrasts] = tfeQCMStimuliToDirectionsContrasts(LMSContrastMat, ...
            'zeroContrastDirection',indDirectionDirections(:,1),'precision',directionPrecision);
        
        % Take the median across voxels
        theMedianTimeCourse = median(fullCleanData(:,:,(jj+((sessionNum-1)*10))),1);
        
        % Apply high pass filter if flag set to true
        if p.Results.highpass
            theMedianTimeCourse = highpass(theMedianTimeCourse,5/288,1/.8);
        end
        
        %%  Make the IAMP packet
        thePacket.response.values   = theMedianTimeCourse;
        thePacket.response.timebase = stimulusStruct.timebase;
        % the stimulus
        thePacket.stimulus.timebase = stimulusStruct.timebase;
        thePacket.stimulus.values   = stimulusStruct.values;
        % the kernel
        thePacket.kernel = analysisParams.HRF;
        
        % Apply high pass filter to HRF if flag set to true
        if p.Results.highpass
            thePacket.kernel.values = highpass(thePacket.kernel.values ,5/288,1/.8);
        end
        
        % the meta data (this is the constrast and directions)
        thePacket.metaData.stimDirections = stimDirections;
        thePacket.metaData.stimContrasts  = stimContrasts;
        thePacket.metaData.lmsContrast    = LMSContrastMat;
        
        count = ((sessionNum-1)*analysisParams.numAcquisitions) + jj;
        iampTimeCoursePacketPocket{count} = thePacket;
        
    end
   
end

end


