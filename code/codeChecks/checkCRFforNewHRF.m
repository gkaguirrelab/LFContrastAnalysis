% this is a script to plot the beta weights that come out of geoffs gear as a CRF

%% This is a demo to help us fit the QCM directly to the time course
%

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('KAS25');

% Set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

% Turn on or off plotting
analysisParams.showPlots = true;

% Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;

% Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

% set the HRF
% load(fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'),'LFContrastAnalysis','subjectHRFs',analysisParams.expSubjID,[analysisParams.expSubjID '_eventGain_results.mat']));
% xBase = zeros(1,analysisParams.expLengthTR);
% xBase(1:length(results.hrf')) = results.hrf';
% analysisParams.HRF.values = xBase;
% analysisParams.HRF.timebase =   analysisParams.timebase*1000;

analysisParams.HRF = generateHRFKernel(6,12,10,analysisParams.timebase*1000);

%% CURRENT WAY -  get the data and the beta weight from my analysis code
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','averageBaseline');
end

averageIampParams = iampOBJ.averageParams(concatParams);

%% USING THE GEAR
% get the beta weights from the results of Geoff's gear
gearParams = results.params(results.meta.vxs(1),1:41);

scale = averageIampParams.paramMainMatrix'/gearParams

gearIampParams.paramNameCell   = {'amplitude'};
gearIampParams.paramMainMatrix = gearParams'*scale;
gearIampParams.matrixRows      = 41;
gearIampParams.matrixCols      = 1;
gearIampParams.noiseSd         = 0;

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,gearIampParams);


% Fit the CRF -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket});

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});

% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict CRF from direction model fits

% Predict the responses for CRF with params from NR common Amp.
crfPlot.respNrCrf = nrCrfOBJ.computeResponse(nrCrfParams{1},crfStimulus,[]);
crfPlot.respNrCrf.color = [.5, .3, .8];

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmCrfMeanOBJ.computeResponse(qcmCrfMeanParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.color = [0, 1, 0];

iampOBJ = tfeIAMP('verbosity','none');
[iampPoints, iampSEM] = iampOBJ.averageParams({gearIampParams});
crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampPoints,iampSEM,'subtractBaseline', true);
