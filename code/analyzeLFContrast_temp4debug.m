%% Set up params
% Set the subject
subjId = 'KAS25'; 

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);

% Set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

% Turn on or off plotting
analysisParams.showPlots = true;

%% Load the relevant data (SDM, HRF, TC)
% Get stimulus design matrix for the entire measurment set (session 1 and session 2 pair)
[stimCells] = makeStimMatrices(subjId); 

%set the HRF
[analysisParams] = loadHRF(analysisParams);

% Load the time course
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

%% Concatenate the data
% Pull out the median time courses
[analysisParams, iampTimeCoursePacketPocket, ~, ~, ~, rawTC] = fit_IAMP(analysisParams,fullCleanData,'concatAndFit', true);

% Concat the stim matrices and time courses
theSignal = [rawTC{1}.values, rawTC{2}.values];
theStimIAMP   =  cat(2, stimCells{:});

% Create timebase
numTimePoints = length(theSignal);
timebase = linspace(0,(numTimePoints-1)*analysisParams.TR,numTimePoints)*1000;

%% Create the IAMP packet for the full experiment session 1 and 2
% Tull time course
thePacketIAMP.response.values   = theSignal;

% Timebase 
thePacketIAMP.response.timebase = timebase;

% Stimulus design matrix remove attentional events)
thePacketIAMP.stimulus.values   = theStimIAMP(1:end-1,:);

% Timebase
thePacketIAMP.stimulus.timebase = timebase;

% The Kernel
kernelVec = zeros(size(timebase));
kernelVec(1:length(analysisParams.HRF.values)) = analysisParams.HRF.values;
thePacketIAMP.kernel.values = kernelVec;
thePacketIAMP.kernel.timebase = timebase;

% Packet meta data
thePacketIAMP.metaData = [];

% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(thePacketIAMP.stimulus.values,1);
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(thePacketIAMP,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

%% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
theStimQCM   =  [directionTimeCoursePacketPocket{1}.stimulus.values,directionTimeCoursePacketPocket{2}.stimulus.values];

% Create the time course packet
timeCoursePacket.response = thePacketIAMP.response;
timeCoursePacket.stimulus.values   = theStimQCM;
timeCoursePacket.stimulus.timebase = timebase;
% the kernel
timeCoursePacket.kernel = thePacketIAMP.kernel;
% packet meta data
timeCoursePacket.metaData = [];

% ###### FIX ###################
% remove subraction of the baseline
% ##############################

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampParams);
% Fit the direction based models to the mean IAMP beta weights

% Fit the time course -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {timeCoursePacket});

% Fit the CRF
[qcmCrfOBJ,qcmCrfParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket},'fitErrorScalar',1000);

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmTimeCourseOBJ,qcmTimeCourseParams] = fitDirectionModel(analysisParams, 'qcmFit', {timeCoursePacket},'fitErrorScalar',1000);

% Fit the CRF with the QCM -- { } is because this expects a cell
% --INITIALIZED
[~,qcmTimeCourseInitParams] = fitDirectionModel(analysisParams, 'qcmFit', {timeCoursePacket},'fitErrorScalar',1000, 'intialParams', qcmCrfParams);

% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict CRF from direction model fits

% Predict the responses for CRF with params from NR common Amp.
crfPlot.respNrCrf = nrCrfOBJ.computeResponse(nrCrfParams{1},crfStimulus,[]);
crfPlot.respNrCrf.color = [.5, .3, .8];

% Predict the responses for CRF with params from QCM fit to CRF stim 
crfPlot.respQCMCrf = qcmTimeCourseOBJ.computeResponse(qcmCrfParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.color = [0, 1, 0];

% Predict the responses for CRF with params from QCM fit to time course
crfPlot.respQCMfullTC = qcmTimeCourseOBJ.computeResponse(qcmTimeCourseParams{1},crfStimulus,[]);
crfPlot.respQCMfullTC.color = [.1, 0, .5];

% Predict the responses for CRF with params from QCM fit to time course
crfPlot.respQCMfullTCInit = qcmTimeCourseOBJ.computeResponse(qcmTimeCourseInitParams{1},crfStimulus,[]);
crfPlot.respQCMfullTCInit.color = [.2, .2, .2];

% dummy up sem as 0s
% THIS NEEDS TO BE CALCULATED WITH THE BOOTSTRAP ROUTINE BUT FOR NOW ZERO
% TO HAVE SOMETHING TO PLOT AND NOT CRASH
iampParams.paramMainMatrix(end) = [];
iampParams.matrixRows = size(iampParams.paramMainMatrix,1);
semParams = iampParams;
semParams.paramMainMatrix =zeros(size(iampParams.paramMainMatrix));

% Plot the CRF from the IAMP, QCM, and  fits
crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semParams,'subtractBaseline', true);


% Get the time course predicitions of the CRF params

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmp = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParams{1}, {timeCoursePacket}, 'plotColor', [0, 0, 1]);

% Get the time course predicitions fromt the QCM params fit to theCRF
timeCoursePlot.qcmCRF = responseFromPacket('qcmPred', analysisParams, qcmCrfParams{1}, {timeCoursePacket}, 'plotColor', [1, 0, 0]);

% Get the time course predicitions fromt the QCM params fit to the time coruse
timeCoursePlot.qcmFullTC = responseFromPacket('qcmPred', analysisParams, qcmTimeCourseParams{1}, {timeCoursePacket}, 'plotColor', [0, 1, 0]);


% Add clean time
timeCoursePlot.timecourse = {thePacketIAMP.response};
timeCoursePlot.timecourse{1}.plotColor =[0, 0, 0];

% Plot the time course prediction for each run using the different fits to
% the crf
plotTimeCourse(analysisParams, timeCoursePlot, 0, 1);



