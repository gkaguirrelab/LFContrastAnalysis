% Set the subject
subjId = 'KAS25'; 

%Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);

%set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

%turn on or off plotting
analysisParams.showPlots = true;

%Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;
%Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

%% Get stimulus design matrix for the entire measurment set (session 1 and session 2 pair)
[stimCells] = makeStimMatrices(subjId); 

%set the HRF
[analysisParams] = loadHRF(analysisParams);

% Load the time course
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Pull out the median time courses
[analysisParams, iampTimeCoursePacketPocket, ~, ~, ~, rawTC] = fit_IAMP(analysisParams,fullCleanData,'concatAndFit', true);

% Concat the stim matrices and time courses
theSignal = [rawTC{1}.values, rawTC{2}.values];
theStimIAMP   =  cat(2, stimCells{:});

% Create timebase
numTimePoints = length(theSignal);
timebase = linspace(0,(numTimePoints-1)*analysisParams.TR,numTimePoints)*1000;

% Create the IAMP packet for the full experiment session 1 and 2
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
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(thePacketIAMP,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

% Get directon/contrast form of time course and IAMP crf packet pockets.
% 
% This conversion is possible because the IAMP packet pocket has meta data
% that we put there to allow exactly this conversion.  That meta data
% encapsulates the key things we need to know about the stimulus obtained
% from the analysis parameters.
% directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
% directionTimeCoursePacketPocket = {directionTimeCoursePacketPocket{1,:},directionTimeCoursePacketPocket{2,:}};

%% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
theStimQCM   =  [directionTimeCoursePacketPocket{1}.stimulus.values,directionTimeCoursePacketPocket{2}.stimulus.values];

% Create the packet
thePacket.response = thePacketIAMP.response;

thePacket.stimulus.values   = theStimQCM;
thePacket.stimulus.timebase = timebase;
% the kernel
thePacket.kernel = thePacketIAMP.kernel;
% packet meta data
thePacket.metaData = [];

% This puts together pairs of acquistions from the two sessions, so that
% we have one IAMP fit for each pair.  We do this because to fit the
% quadratic model, we need data for all of the color directions together.
% 
% NOTE: This bit is very specific to the design of the experiment we are
% currently analyzing, and has to do specifically with the way color
% directions were studied across acquisitions and sessions.
% 
% 
% ###### FIX ###################
% remove subraction of the baseline
% ##############################

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampParams);
%remove attentional event regressor
directionCrfMeanPacket.response.values(end) = []; 
% Fit the direction based models to the mean IAMP beta weights

% Fit the CRF -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {thePacket});

% % Fit the CRF with the NR common amplitude -- { } is because this expects a cell
% [~,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', {thePacket}, 'commonAmp', true);
% 
% % Fit the CRF with the NR common Exponent -- { } iPs because this expects a cell
% [~,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', {thePacket}, 'commonExp', true);
% 
% % Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
% [~,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', {thePacket}, 'commonAmp', true, 'commonExp', true);

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {thePacket},'fitErrorScalar',1000);

% Do some plotting of these fits
if analysisParams.showPlots
    [nrVals] = plotNakaRushtonFromParams(qcmCrfMeanParams{1}.crfAmp ,qcmCrfMeanParams{1}.crfExponent,qcmCrfMeanParams{1}.crfSemi,...
        'analysisParams',analysisParams,'plotFunction',true,'savePlot',true);
end
% ######### FIX ###################
% make a function that subtracts the baseline for generating the crf plots
% ################################


% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict CRF from direction model fits

% Predict the responses for CRF with params from NR common Amp.
crfPlot.respNrCrf = nrCrfOBJ.computeResponse(nrCrfParams{1},crfStimulus,[]);
crfPlot.respNrCrf.color = [.5, .3, .8];

% % Predict the responses for CRF with params from NR common Amp.
% crfPlot.respNrCrfAmp = nrCrfOBJ.computeResponse(nrCrfParamsAmp{1},crfStimulus,[]);
% crfPlot.respNrCrfAmp.color = [0, 0, 1];
% 
% % Predict the responses for CRF with params from NR common Exp
% crfPlot.respNrCrfExp = nrCrfOBJ.computeResponse(nrCrfParamsExp{1},crfStimulus,[]);
% crfPlot.respNrCrfExp.color = [0, .33, 1];
% 
% % Predict the responses for CRF with params from NR common Amp and Exps
% crfPlot.respNrCrfAmpExp = nrCrfOBJ.computeResponse(nrCrfParamsAmpExp{1},crfStimulus,[]);
% crfPlot.respNrCrfAmpExp.color = [0, .66, 1];

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmCrfMeanOBJ.computeResponse(qcmCrfMeanParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.color = [0, 1, 0];

% % Now use the QCM to get NR parameters that can be applied to crfStimulus using the
% % Naka-Rushton objects.
% nrDirections = nrCrfOBJ.directions;
% nrContrasts = ones(1,size(nrDirections,2));
% tempStimulus.values = [nrDirections ; nrContrasts];
% tempStimulus.timebase = 1:size(nrDirections,2);
% tempResp = qcmCrfMeanOBJ.computeResponse(qcmCrfMeanParams{1},tempStimulus,[]);
% for ii = 1:length(nrCrfParamsAmpExp{1})
%     nrQcmBasedParams{1}(ii) = nrCrfParamsAmpExp{1}(ii);
%     nrQcmBasedParams{1}(ii).crfSemi = qcmCrfMeanParams{1}.crfSemi/tempResp.metaData.quadraticFactors(ii);
%     nrQcmBasedParams{1}(ii).crfExponent = qcmCrfMeanParams{1}.crfExponent;
%     nrQcmBasedParams{1}(ii).crfAmp = qcmCrfMeanParams{1}.crfAmp;
%     nrQcmBasedParams{1}(ii).crfOffset = qcmCrfMeanParams{1}.crfOffset;
% end
% crfPlot.respNrQcmBased = nrCrfOBJ.computeResponse(nrQcmBasedParams{1},crfStimulus,[]);
% crfPlot.respNrQcmBased.color = [1 0 0];
% 
% % Fit the CRF with the NR common amplitude and semisaturation  -- { } is because this expects a cell
% % This time start with parameters unpacked from QCM filt.
% [nrQcmBasedCrfOBJ,nrQcmBasedCrfParamsAmpSemi] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, ...
%     'commonAmp', true, 'commonExp', true, 'initialParams', nrQcmBasedParams{1});
% crfPlot.respNrQcmBasedCrfAmpSemi = nrCrfOBJ.computeResponse(nrQcmBasedCrfParamsAmpSemi{1},crfStimulus,[]);
% crfPlot.respNrQcmBasedCrfAmpSemi.color = [1 0.2 0];

% dummy up sem as 0s
% THIS NEEDS TO BE CALCULATED WITH THE BOOTSTRAP ROUTINE BUT FOR NOW RAND
% TO HAVE SOMETHING TO PLOT AND NOT CRASH
iampParams.paramMainMatrix(end) = [];
iampParams.matrixRows = size(iampParams.paramMainMatrix,1)
semParams = iampParams;
semParams.paramMainMatrix =zeros(size(iampParams.paramMainMatrix));

% Plot the CRF from the IAMP, QCM, and  fits
if analysisParams.showPlots
    %[iampPoints, iampSEM] = iampOBJ.averageParams(concatParams);
    crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semParams,'subtractBaseline', true);
    figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_CRF_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    FigureSave(figNameCrf,crfHndl,'pdf');
end

% Get the time course predicitions of the CRF params

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmp = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmp{1}, directionTimeCoursePacketPocket, 'plotColor', [0, 0, 1]);

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmpSemi = responseFromPacket('nrPred', analysisParams, nrCrfParamsExp{1}, directionTimeCoursePacketPocket, 'plotColor', [0, .33, 1]);

% Get the time course predicitions from the NR common Amp and Semi fit to the CRF
timeCoursePlot.nrAmpSemiExp = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmpExp{1}, directionTimeCoursePacketPocket, 'plotColor', [0, .66, 1]);

% Get the time course predicitions fromt the QCM params fit to the CRF
timeCoursePlot.qcm = responseFromPacket('qcmPred', analysisParams, qcmCrfMeanParams{1}, directionTimeCoursePacketPocket, 'plotColor', [0, 1, 0]);

% Get the time course predicitions from the NR common Amp and Semi fit to
% the CRF, based on QCM fit.
timeCoursePlot.nrQcmBasedAmpSemi = responseFromPacket('nrPred', analysisParams, nrQcmBasedCrfParamsAmpSemi{1}, directionTimeCoursePacketPocket, 'plotColor', [0.5 0.2 0.6]);

% Get the predictions from individual IAMP params
for ii = 1:size(iampParams,1)
    for jj = 1:size(iampParams,2)
        timeCoursePlot.iamp{ii,jj} = responseFromPacket('IAMP', analysisParams, iampParams, iampTimeCoursePacketPocket{ii,jj}, 'plotColor', [0.5 0.2 0]);
    end
end

% FIX THIS
% Get the time course prediction from the avarage IAMP params
iampParamsTC.sessionOne  = iampOBJ.averageParams(iampParams(1,:));
iampParamsTC.sessionTwo  = iampOBJ.averageParams(iampParams(2,:));
iampParamsTC.baseline = concatBaselineShift;
timeCoursePlot.meanIamp = responseFromPacket('meanIAMP', analysisParams, iampParamsTC, iampTimeCoursePacketPocket, 'plotColor', [0.0 0.2 0.6]);

% Add clean time
timeCoursePlot.timecourse = rawTC;


% Plot the time course prediction for each run using the different fits to
% the crf
if analysisParams.showPlots
    tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, concatBaselineShift, analysisParams.numSessions*analysisParams.numAcquisitions);
    figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    FigureSave(figNameTc,tcHndl,'pdf');
end

% Plot isoresponce contour
if analysisParams.showPlots
    thresholds = [0.10, 0.2, 0.3];
    colors     = [0.5,0.0,0.0; 0.5,0.5,0.0; 0.0,0.5,0.5;];
    qcmHndl    = plotIsoresponse(analysisParams,iampPoints,qcmCrfMeanParams,thresholds,nrCrfParamsAmp,colors);
    figNameQcm = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_QCM_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    FigureSave(figNameQcm,qcmHndl,'pdf');
end
