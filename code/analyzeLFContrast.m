%% Set up params
% Set the subject: 'LZ23', 'KAS25', 'AP26'
subjId = 'KAS25'; 

% Load the subject relevant info
analysisParams = getSubjectParams(subjId);

%set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

%turn on or off plotting
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
% Full time course
thePacketIAMP.response.values   = theSignal;

% Timebase
thePacketIAMP.response.timebase = timebase;

% Stimulus design matrix remove attentional events
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

% Create the tine course packet
thePacket.response = thePacketIAMP.response;
thePacket.stimulus.values   = theStimQCM;
thePacket.stimulus.timebase = timebase;
thePacket.kernel = thePacketIAMP.kernel;
thePacket.metaData = [];

% Make the CRF packet
directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampParams);

%% Fit the direction based models to the mean IAMP beta weights
% Fit the CRF -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {thePacket});

% Fit the CRF with the NR common amplitude -- { } is because this expects a cell
[~,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', {thePacket}, 'commonAmp', true);

% Fit the CRF with the NR common Exponent -- { } iPs because this expects a cell
[~,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', {thePacket}, 'commonExp', true);

% Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
[~,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', {thePacket}, 'commonAmp', true, 'commonExp', true);

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {thePacket},'fitErrorScalar',1000);

% Do some plotting of these fits
if analysisParams.showPlots
    [nrVals] = plotNakaRushtonFromParams(qcmCrfMeanParams{1}.crfAmp ,qcmCrfMeanParams{1}.crfExponent,qcmCrfMeanParams{1}.crfSemi,...
        'analysisParams',analysisParams,'plotFunction',true,'savePlot',true);
end

% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict CRF from direction model fits

% Predict the responses for CRF with params from NR common Amp.
crfPlot.respNrCrf = nrCrfOBJ.computeResponse(nrCrfParams{1},crfStimulus,[]);
crfPlot.respNrCrf.color = [.5, .3, .8];

% Predict the responses for CRF with params from NR common Amp.
crfPlot.respNrCrfAmp = nrCrfOBJ.computeResponse(nrCrfParamsAmp{1},crfStimulus,[]);
crfPlot.respNrCrfAmp.color = [0, 0, 1];

% Predict the responses for CRF with params from NR common Exp
crfPlot.respNrCrfExp = nrCrfOBJ.computeResponse(nrCrfParamsExp{1},crfStimulus,[]);
crfPlot.respNrCrfExp.color = [0, .33, 1];

% Predict the responses for CRF with params from NR common Amp and Exps
crfPlot.respNrCrfAmpExp = nrCrfOBJ.computeResponse(nrCrfParamsAmpExp{1},crfStimulus,[]);
crfPlot.respNrCrfAmpExp.color = [0, .66, 1];

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmCrfMeanOBJ.computeResponse(qcmCrfMeanParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.color = [0, 1, 0];

% dummy up sem as 0s
% THIS NEEDS TO BE CALCULATED WITH THE BOOTSTRAP ROUTINE BUT FOR NOW RAND
% TO HAVE SOMETHING TO PLOT AND NOT CRASH
iampParams.paramMainMatrix(end) = [];
iampParams.matrixRows = size(iampParams.paramMainMatrix,1)
semParams = iampParams;
semParams.paramMainMatrix =zeros(size(iampParams.paramMainMatrix));

% Plot the CRF from the IAMP, QCM, and  fits
if analysisParams.showPlots
    crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semParams,'subtractBaseline', true);
%     figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
%         [analysisParams.expSubjID,'_CRF_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
%     FigureSave(figNameCrf,crfHndl,'pdf');
end

% Get the time course predicitions of the CRF params

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmp = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParamsAmp{1}, {timeCoursePacket}, 'plotColor', [0, 0, 1]);

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmpSemi = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParamsExp{1}, {timeCoursePacket}, 'plotColor', [0, .33, 1]);

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmpSemiExp = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParamsAmpExp{1}, {timeCoursePacket}, 'plotColor', [0, .66, 1]);

% Get the time course predicitions fromt the QCM params fit to the CRF
timeCoursePlot.qcm = responseFromPacket('qcmPred', analysisParams, qcmCrfMeanParams{1}, {timeCoursePacket}, 'plotColor', [0, 1, 0]);

% Add clean time
timeCoursePlot.timecourse = {thePacketIAMP.response};
timeCoursePlot.timecourse{1}.plotColor =[0, 0, 0];

% Plot the time course prediction for each run using the different fits to
% the crf
if analysisParams.showPlots
    tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, 0, 1);
%     figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
%         [analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
%     FigureSave(figNameTc,tcHndl,'pdf');
end

% Plot isoresponce contour
if analysisParams.showPlots
    thresholds = [0.15];
    colors     = [0.5,0.0,0.0];
    qcmHndl    = plotIsoresponse(analysisParams,iampParams,qcmCrfMeanParams,thresholds,nrCrfParamsAmp,colors);
%     figNameQcm = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
%         [analysisParams.expSubjID,'_QCM_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
%     FigureSave(figNameQcm,qcmHndl,'pdf');
end
