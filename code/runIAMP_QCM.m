function [analysisParams,paramsQCMFit, meanIAMPBetas, semIAMPBetas,packetPocket,paramsFitIAMP, fitResponseStructQCM] = runIAMP_QCM(analysisParams,fullCleanData)
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


    %% Construct the model object
    temporalFit = tfeIAMP('verbosity','none');

    %% Create a cell of stimulusStruct (one struct per run)
    for jj = 1:analysisParams.numAcquisitions

        % identify the data param file
        dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});

        % We are about to load the data param file. First silence the warning
        % for EnumberableClass. Save the warning state.
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

        % make stimulus values
        % Stim coding: 80% = 1, 40% = 2, 20% = 3, 10% = 4, 5% = 5, 0% = 6;
        stimulusStruct.values =  createRegressors(expParams,analysisParams.baselineCondNum,totalTime,deltaT);

        %[ * NOTE: MB: make sure the timestep is loaded from the pulse params
        %istead of set here]
        responseStruct.timeStep = analysisParams.timeStep;

        % get attention event regressor
        [~, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', stimulusStruct.timebase);

        % add attention events to regressor matrix
        stimulusStruct.values = [stimulusStruct.values;eventsRegressor];

        % Set the number of instances.
        defaultParamsInfo.nInstances = size(stimulusStruct.values,1);

        % Get kernel
        kernelStruct = generateHRFKernel(6,12,10,stimulusStruct.timebase);

        % make the stimulus portion of packet for fitting
        thePacket.stimulus.timebase = stimulusStruct.timebase;
        thePacket.stimulus.values   = stimulusStruct.values;

        % add the response field
        thePacket.response.timebase = stimulusStruct.timebase;
        timeCourseValues(:,jj,sessionNum) = median(fullCleanData(:,:,(jj+((sessionNum-1)*10))),1)';
        thePacket.response.values =timeCourseValues(:,jj,sessionNum)';
        % add the kernel field
        thePacket.kernel = kernelStruct;

        % add a metaData field
        thePacket.metaData = [];

        %% Perform the fit
        [paramsFit,fVal,IAMPResponses] = ...
            temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, ...
            'searchMethod','linearRegression');

        % Plot data and IAMP fit
        if(analysisParams.generateIAMPPlots)
            temporalFit.plot(thePacket.response,'Color',[1 0 0]);
            temporalFit.plot(IAMPResponses,'Color',[0 1 0],'NewWindow',false);
        end
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
end
meanIAMPBetas = [meanIAMPBetas;meanBaseline];
semIAMPBetas  = [semIAMPBetas;semBaseline];

% ADD CROSS VALIDATION HERE
[rmseMeanIAMP rmseMeanQCM rmseSemIAMP rmseSemQCM] = crossValidateIAMP_QCM(analysisParams,betas,timeCourseValues)

%% Fit IAMP crfs with QCM
% Set parameters and construct a QCM object.
temporalFitQCM = tfeQCM('verbosity','none','dimension',analysisParams.theDimension);

%% Set up contrast values matched to resoponse order
% Set up stim order info to creat LMS contrast by timepoint matrix
stimulusStruct.values   = [generateStimCombinations(analysisParams.contrastCoding,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension),[0;0]];
stimulusStruct.timebase = 1:length(stimulusStruct.values);

%% Snag response values from IAMP fit.
thePacket.response.values = meanIAMPBetas';
thePacket.response.timebase = 1:length(thePacket.response.values);

%% Construct a packet for the QCM to fit.
thePacket.stimulus = stimulusStruct;
thePacket.kernel = [];
thePacket.metaData = [];

%% Fit
% allow QCM to fit the offset 
defaultParamsInfo.noOffset = false;
[paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket,'defaultParamsInfo',defaultParamsInfo);
fprintf('Model parameter from fits:\n');
temporalFitQCM.paramPrint(paramsQCMFit)

end
