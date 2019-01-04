function [analysisParams,paramsQCMFit, meanIAMPBetas, semIAMPBetas,packetPocket,paramsFitIAMP, fitResponseStructQCM, fitParams, meanNRParams, crossval] = runIAMP_QCM(analysisParams,fullCleanData)
% Takes in the clean time series data and the analysis params and runs the IMAP-QCM model.
%
% Syntax:
%   [analysisParams,paramsQCMFit, meanIAMPBetas, semIAMPBetas,packetPocket, ...
%    paramsFitIAMP, fitResponseStructQCM, baselineBetas] = runIAMP_QCM(analysisParams,fullCleanData)
%
% Description:
%    This function takes in the clean time series data and the analysis params
%    and runs the IMAP-QCM model. This function builds a stimulus design matirx
%    based on the analysisParams (from each run of the experiemnt) and run the
%    IAMP model on the cleaned and trial sorted data. The beta weights from the
%    IAMP model are then fit with the QCM model.
%
%
% Inputs:
%    inFile            - File name of a text file. (string)
%
% Outputs:
%    analysisParams          - Returns analysisParams with any updates
%    paramsQCMFit            - Paramerter fits of the QCM to the IAMP betas
%    meanIAMPBetas           - Mean of the IAMP beta weights across all runs.
%    semIAMPBetas            - Stanard error of mean beta weights
%    packetPocket            - Struct with the IAMP packet for each run
%    paramsFitIAMP           - Parameter fits for the IAMP model
%    fitResponseStructQCM    - QCM Fit struct
%    baselineBetas           - beta weigths for the baseline condition fro each run
%
% Optional key/value pairs:
%    none

% MAB 09/09/18

IAMPBetas = [];
IAMPsem  = [];
baselineBetas = [];
count = 1;
for sessionNum = 1:length(analysisParams.sessionFolderName)
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'melaDataPath'), analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderFiles = textFile2cell(trialOrderFile);
    
    % Get the Directions of each session. This requires analysisParams.directionCoding to be organized
    % such that the directions are grouped per session and these groups are in the same  order as the
    % sessions order
    
    % get then number of direction over all sessiosn
    analysisParams.numDirections = size(unique(analysisParams.directionCoding','rows')',2);
    
    % check how the sessions were ordered.
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
    temporalFit = tfeIAMP('verbosity','none');
    
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
        indDirectionDirections = round(directionCoding(1:2,:),directionPrecision);
        LMSContrastMat(3,:) = [];
        [stimDirections,stimContrasts] = tfeQCMStimuliToDirectionsContrasts(LMSContrastMat, ...
            'zeroContrastDirection',indDirectionDirections(:,1),'precision',directionPrecision);
        
        % Get the unique directions and remove the [0;0] column
        checkDirections = tfeQCMParseDirections(stimDirections);
        % [WRITE A CHECK HERE THAT EACH COLUMN OF indDirectionDirections is
        % in in checkDirections, and vice versa.
        
        %[ * NOTE: MB: make sure the timestep is loaded from the pulse params
        %istead of set here]
        responseStruct.timeStep = analysisParams.timeStep;
        
        % get attention event regressor
        [~, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', stimulusStruct.timebase);
        
        % add attention events to regressor matrix
        stimulusStruct.values = [stimulusStruct.values;eventsRegressor];
        
        % Set the number of instances.
        clear defaultParamsInfo
        defaultParamsInfo.nInstances = size(stimulusStruct.values,1);
        
        % Get kernel
        kernelStruct = generateHRFKernel(6,12,10,stimulusStruct.timebase);
        
        % make the stimulus portion of packet for fitting
        thePacket.stimulus.timebase = stimulusStruct.timebase;
        thePacket.stimulus.values   = stimulusStruct.values;
        
        % add the response field
        responses = median(fullCleanData(:,:,(jj+((sessionNum-1)*10))),1)';
        thePacket.response.timebase = stimulusStruct.timebase;
        timeCourseValues(:,jj,sessionNum) = responses;
        thePacket.response.values =timeCourseValues(:,jj,sessionNum)';
        % add the kernel field
        thePacket.kernel = kernelStruct;
        
        % add a metaData field
        thePacket.metaData = [];
        
        % Perform the fit
        [paramsFit,fVal,IAMPResponses] = ...
            temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, ...
            'searchMethod','linearRegression');
        
        % Look at IAMP fit contrast response functions
        %         stimulusStruct.values   = [generateStimCombinations(analysisParams.contrastCoding,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension),[0;0]];
        %         stimulusStruct.timebase = 1:length(stimulusStruct.values);
        %         figure; clf; hold on;
        %         plot(paramsFit.paramMainMatrix(1:end-1));
        
        % Fit Naka-Ruston to the timecourse with things common across directions
        NOOFFSET = false;
        commonAmp = false;
        commonSemi = false;
        commonExp = false;
        commonOffset = true;
        
        % Create the NR Pactket
        theNRPacket.stimulus.values   = [stimDirections ; stimContrasts];
        theNRPacket.stimulus.timebase = thePacket.stimulus.timebase;
        theNRPacket.response          = thePacket.response;
        theNRPacket.kernel            = thePacket.kernel;
        theNRPacket.metaData          = [];
        
        % Init the tfeNakaRushtonDirection object
        NRDirectionObj = tfeNakaRushtonDirection(indDirectionDirections, ...
            'lockOffsetToZero',NOOFFSET,'commonAmp',commonAmp,'commonSemi',commonSemi,'commonExp',commonExp,'commonOffset',commonOffset);
        
        % Fit the packet
        %[fitNRDirectionParams{sessionNum,jj},~,NRDirectionFitResponses] = NRDirectionObj.fitResponse(theNRPacket);
        %fprintf('\nNaka-Rushton parameters from fit:\n');
        %NRDirectionObj.paramPrint(fitNRDirectionParams{sessionNum,jj});
        
        % Plot data and IAMP fit
        %         if(analysisParams.generateIAMPPlots)
        %             temporalFit.plot(thePacket.response,'Color',[1 0 0]);
        %             temporalFit.plot(IAMPResponses,'Color',[0 1 0],'NewWindow',false);
        %             NRDirectionObj.plot(NRDirectionFitResponses,'Color',[0 0 1],'NewWindow',false);
        %         end
        paramsFitIAMP{count} = paramsFit;
        packetPocket{count} = thePacket;
        
        % Remove the meta weight for the attentional event
        betas(:,jj,sessionNum)= paramsFit.paramMainMatrix(1:end-1);
        count = count+1;
    end
    
    % Calculate mean of the betas
    IAMPBetas = [IAMPBetas, mean(betas(1:end-1,:,sessionNum),2)];
    IAMPsem = [IAMPsem, std(betas(1:end-1,:,sessionNum),0,2)./sqrt(analysisParams.numAcquisitions)];
end

% Get the mean Naka-Ruston paramters for each direction (as fit with the
% TFE NR)
% for kk = 1:size(fitNRDirectionParams,1)
%     for pp = 1:size(fitNRDirectionParams,2)
%
%         tmpParams = fitNRDirectionParams{kk,pp};
%
%         for gg = 1:size(tmpParams,2);
%             crfAmp(pp,gg)      = tmpParams(gg).crfAmp;
%             crfExponent(pp,gg) = tmpParams(gg).crfExponent;
%             crfOffset(pp,gg)   = tmpParams(gg).crfOffset;
%             crfSemi(pp,gg)     = tmpParams(gg).crfSemi;
%             expFalloff(pp,gg)  = tmpParams(gg).expFalloff;
%             noiseSd(pp,gg)     = tmpParams(gg).noiseSd;
%         end
%     end
%
%     meanNRParams.crfAmp(kk,:) = mean(crfAmp);
%     meanNRParams.crfExponent(kk,:) = mean(crfExponent);
%     meanNRParams.crfOffset(kk,:) = mean(crfOffset);
%     meanNRParams.crfSemi(kk,:) = mean(crfSemi);
%     meanNRParams.expFalloff(kk,:) = mean(expFalloff);
%     meanNRParams.noiseSd(kk,:) = mean(noiseSd);
%
% end

baselineBetas = betas(end,:,:);
meanBaseline = mean(baselineBetas(:));
semBaseline = std(baselineBetas(:))./sqrt(numel(betas(end,:,:)));

% combine the betas across sessions
meanIAMPBetas = [];
semIAMPBetas  = [];
baselineCond  = [];
for pp = 1:size(IAMPBetas,2)
    meanIAMPBetas = [meanIAMPBetas; IAMPBetas(:,pp)];
    semIAMPBetas  = [semIAMPBetas; IAMPsem(:,pp)];
    fitParams.meanIAMPBetas{pp} = paramsFit;
    fitParams.meanIAMPBetas{pp}.paramMainMatrix = [];
    fitParams.meanIAMPBetas{pp}.paramMainMatrix = IAMPBetas(:,pp);   
end
meanIAMPBetas = [meanIAMPBetas;meanBaseline];
semIAMPBetas  = [semIAMPBetas;semBaseline];


% ADD CROSS VALIDATION HERE
% [rmseMeanIAMP, rmseMeanQCM, rmseSemIAMP, rmseSemQCM] = crossValidateIAMP_QCM(analysisParams,betas,timeCourseValues, paramsFitIAMP, packetPocket, 'showPlots', true);
% crossval = [rmseMeanIAMP rmseMeanQCM rmseSemIAMP rmseSemQCM];

end


