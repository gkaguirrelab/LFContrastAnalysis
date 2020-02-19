function [modelResponseStructIAMP, modelResponseStructQCM, paramsFitIAMP, paramsFitQCM] = fullTimeCourseFit(subjId)


% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);

% set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

% turn on or off plotting
analysisParams.showPlots = true;

% Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;
% Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

%set the HRF
[analysisParams] = loadHRF(analysisParams);


%% Get stimulus design matrix for the entire measurment set (session 1 and session 2 pair)
[stimCells] = makeStimMatrices(subjId);

% Get the time course data
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Pull out the median time courses
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData,'concatAndFit', true);

% Concat the stim matrices and time courses
theSignal = [rawTC{1}.values, rawTC{2}.values];
theStimIAMP   =  cat(2, stimCells{:});

% Create timebase
numTimePoints = length(theSignal);
timebase = linspace(0,(numTimePoints-1)*analysisParams.TR,numTimePoints)*1000;

% Create the packet
thePacketIAMP.response.values   = theSignal;
thePacketIAMP.response.timebase = timebase;

thePacketIAMP.stimulus.values   = theStimIAMP;
thePacketIAMP.stimulus.timebase = timebase;
% the kernel
kernelVec = zeros(size(timebase));
kernelVec(1:length(analysisParams.HRF.values)) = analysisParams.HRF.values;
thePacketIAMP.kernel.values = kernelVec;
thePacketIAMP.kernel.timebase = timebase;
% packet meta data
thePacketIAMP.metaData = [];

% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(thePacketIAMP.stimulus.values,1);
[paramsFitIAMP,fVal,IAMPResponses] = iampOBJ.fitResponse(thePacketIAMP,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
% generate time course from params fit and stim struct
modelResponseStructIAMP = iampOBJ.computeResponse(paramsFitIAMP,thePacketIAMP.stimulus,thePacketIAMP.kernel);

%% Do the QCM fit!
% Get directon/contrast form of time course and IAMP crf packet pockets.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
theStimQCM   =  [directionTimeCoursePacketPocket{1}.stimulus.values,directionTimeCoursePacketPocket{2}.stimulus.values];

% Create the packet
thePacketQCM.response = thePacketIAMP.response;

thePacketQCM.stimulus.values   = theStimQCM;
thePacketQCM.stimulus.timebase = timebase;
% the kernel
thePacketQCM.kernel = thePacketIAMP.kernel;
% packet meta data
thePacketQCM.metaData = [];

% generate time course from params fit and stim struct
[qcmOBJ,paramsFitQCM] = fitDirectionModel(analysisParams, 'qcmFit', {thePacketQCM});
qcmOBJ.paramPrint(paramsFitQCM{1})
modelResponseStructQCM = qcmOBJ.computeResponse(paramsFitQCM{1},thePacketQCM.stimulus,thePacketQCM.kernel);


