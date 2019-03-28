% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26');

% Make mask from the area and eccentricity maps
analysisParams.areaNum     = 1;
analysisParams.eccenRange  = [1 20];

% Define the TR
analysisParams.TR = 0.800;
analysisParams.baselineCondNum = 6;
analysisParams.timeStep = 1/100;
analysisParams.generateIAMPPlots = false;
analysisParams.generateCrossValPlots = false;

% Plotting params
 analysisParams.numSamples = 25;

% Get the cleaned time series
[fullCleanData, analysisParams] = getTimeCourse(analysisParams);

%% Run the IAMP/QCM models
%
% Fit IAMP 
% 
% Fit IAMP to each constructed packet and create packetPocket cell array of
% all the fit packets.
%     packetPocket - Meta data of packePocket contains the direction/contrast form of the same packet.
%     iampOBJ - the tfe IAMP object
%     iampParams - cell array of iampParams for each object
%
% NOTE: Each session gets its own row in the packet pocket.  May want to sweep
% back at some point and match conventions in analysis params to this, for
% example by making the various cell arrays columns rather than rows to
% match.  Similarly with LMVectorAngles vector, which could turn into a
% matrix.
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% Get directon/contrast form of time course and IAMP crf packet pockets.
%
% This conversion is possible because the IAMP packet pocket has meta data
% that we put there to allow exactly this conversion.  That meta data
% encapsulates the key things we need to know about the stimulus obtained
% from the analysis parameters.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);

% This puts together pairs of acquistions from the two sessions, so that
% we have one IAMP fit for each pair.  We do this because to fit the
% quadratic model, we need data for all of the color directions together.
%
% NOTE: This bit is very specific to the design of the experiment we are
% currently analyzing, and has to do specifically with the way color
% directions were studied across acquisitions and sessions.
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampOBJ.averageParams(concatParams));

%% Fit the direction based models to the mean IAMP beta weights 
%
% Fit the CRF with the NR common amplitude -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true);

% Fit the CRF with the NR common Exponent -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonExp', true);

% Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonExp', true);

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});

%% Do some plotting of these fits
%
% Upsample the NR repsonses 
crfStimulus = upsampleCRF(analysisParams);

%% Predict CRF from direction model fits
%
% Predict the responses for CRF with params from NR common Amp.
crfPlot.respNrCrfAmp = nrCrfOBJ.computeResponse(nrCrfParamsAmp{1},crfStimulus,[]);
crfPlot.respNrCrfAmp.color = [0, 0, 1];

% Predict the responses for CRF with params from NR common Exp
crfPlot.respNrCrfExp = nrCrfOBJ.computeResponse(nrCrfParamsExp{1},crfStimulus,[]);
crfPlot.respNrCrfExp.color = [0, .33, 1];

% Predict the responses for CRF with params from NR common Amp and Exp
crfPlot.respNrCrfAmpExp = nrCrfOBJ.computeResponse(nrCrfParamsAmpExp{1},crfStimulus,[]);
crfPlot.respNrCrfAmpExp.color = [0, .66, 1];

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmCrfMeanOBJ.computeResponse(qcmCrfMeanParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.color = [0, 1, 0]; 

%% Now use the QCM to get NR parameters that can be applied to crfStimulus using the
% Naka-Rushton objects.
nrDirections = nrCrfOBJ.directions;
nrContrasts = ones(1,size(nrDirections,2));
tempStimulus.values = [nrDirections ; nrContrasts];
tempStimulus.timebase = 1:size(nrDirections,2);
tempResp = qcmCrfMeanOBJ.computeResponse(qcmCrfMeanParams{1},tempStimulus,[]);
for ii = 1:length(nrCrfParamsAmpExp{1})
    nrQcmBasedParams{1}(ii) = nrCrfParamsAmpExp{1}(ii);
    nrQcmBasedParams{1}(ii).crfSemi = qcmCrfMeanParams{1}.crfSemi/tempResp.metaData.quadraticFactors(ii);
    nrQcmBasedParams{1}(ii).crfExponent = qcmCrfMeanParams{1}.crfExponent;
    nrQcmBasedParams{1}(ii).crfAmp = qcmCrfMeanParams{1}.crfAmp;
    nrQcmBasedParams{1}(ii).crfOffset = qcmCrfMeanParams{1}.crfOffset;
end
% crfPlot.respNrQcmBased = nrCrfOBJ.computeResponse(nrQcmBasedParams{1},crfStimulus,[]);
% crfPlot.respNrQcmBased.color = [1 0 0];

% Fit the CRF with the NR common amplitude and semisaturation  -- { } is because this expects a cell
% This time start with parameters unpacked from QCM filt.
[nrQcmBasedCrfOBJ,nrQcmBasedCrfParamsAmpSemi] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, ...
    'commonAmp', true, 'commonExp', true, 'initialParams', nrQcmBasedParams{1});
crfPlot.respNrQcmBasedCrfAmpSemi = nrCrfOBJ.computeResponse(nrQcmBasedCrfParamsAmpSemi{1},crfStimulus,[]);
crfPlot.respNrQcmBasedCrfAmpSemi.color = [1 0.2 0];

%% Plot the CRF from the IAMP, QCM, and  fits
iampPoints = iampOBJ.averageParams(concatParams);
crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampPoints);
figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
                                [analysisParams.expSubjID,'_CRF_' analysisParams.sessionNickname '.pdf']);
FigureSave(figNameCrf,crfHndl,'pdf');

%% Get the time course predicitions of the CRF params
%
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


% Add clean time 
timeCoursePlot.rawTC = rawTC;


% %Plot the time course prediction for each run using the different fits to
% %the crf

tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, concatBaselineShift, analysisParams.numSessions*analysisParams.numAcquisitions);
figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
                                [analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname '.pdf']);
FigureSave(figNameTc,tcHndl,'pdf');

% % Plot isoresponce contour
thresholds = [0.10, 0.2, 0.3];
colors     = [0.5,0.0,0.0; 0.5,0.5,0.0; 0.0,0.5,0.5;];
qcmHndl    = plotIsoresponse(analysisParams,iampPoints,qcmCrfMeanParams,thresholds,nrCrfParamsAmp,colors);
figNameQcm = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
                                [analysisParams.expSubjID,'_QCM_' analysisParams.sessionNickname '.pdf']);
FigureSave(figNameQcm,qcmHndl,'pdf');
