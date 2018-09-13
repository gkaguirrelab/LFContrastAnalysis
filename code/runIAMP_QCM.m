function [analysisParams,paramsQCMFit, meanIAMPBetas] = runIAMP_QCM(analysisParams,cleanRunData)
% Takes in a text file name and retuns a cell of the lines of the text file
%
% Syntax:
%   filesCell = textFile2cell(inFile)
%
% Description:
%    This function takes in a file name for a trext file and returns a cell
%    that is composed of the lines of the text file. Example of this would
%    be a text file of file names the output is a cell of files names.
%
% Inputs:
%    inFile            - File name of a text file. (string)
%
% Outputs:
%    fileCell          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    none

% MAB 09/09/18

trialOrderDir  = fullfile(getpref(analysisParams.projectName,'melaDataPath'), analysisParams.expSubjID,analysisParams.sessionDate,analysisParams.sessionNumber);
trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName,'experimentFiles','dataFiles.txt');
trialOrderFiles = textFile2cell(trialOrderFile);


%% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

%% Create a cell of stimulusStruct (one struct per run)
for jj = 1:numAcquisitions
    
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
    [~, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', thePacket.response.timebase);
    
    % add attention events to regressor matrix
    stimulusStruct.values = [stimulusStruct.values;eventsRegressor];
    
    % Set the number of instances.
    defaultParamsInfo.nInstances = size(stimulusStruct.values,1);
    
    % Get kernel
    kernelStruct = generateHRFKernel(6,12,10,stimulusStruct.timebase)
    
    % make the stimulus portion of packet for fitting
    thePacket.stimulus.timebase = stimulusStruct.timebase;
    thePacket.stimulus.values   = stimulusStruct.values;
    
    % add the response field
    thePacket.response.timebase = stimulusStruct.timebase;
    thePacket.response.values = median(cleanRunData,1);
    
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
    
    betas(:,jj)= paramsFit.paramMainMatrix;
end

% Calculate mean of the betas
meanIAMPBetas = mean(betas,2);

%% Fit IAMP crfs with QCM
% Set parameters and construct a QCM object.
temporalFitQCM = tfeQCM('verbosity','none','dimension',theDimension);

%% Set up contrast values matched to resoponse order
% Set up stim order info to creat LMS contrast by timepoint matrix
stimulusStruct.values   = [generateStimCombinations(analysisParams.contrastCoding,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension),[0;0]];
stimulusStruct.timebase = 1:length(stimulusStruct.values);

%% Snag response values from IAMP fit.
thePacket.response.values = meanIAMPBetas(1:21)';
thePacket.response.timebase = 1:length(thePacket.response.values);

%% Construct a packet for the QCM to fit.
thePacket.stimulus = stimulusStruct;
thePacket.kernel = [];
thePacket.metaData = [];

%% Fit
[paramsQCMFit,fVal,fitResponseStruct] = temporalFitQCM.fitResponse(thePacket);
fprintf('Model parameter from fits:\n');
temporalFitQCM.paramPrint(paramsQCMFit)



