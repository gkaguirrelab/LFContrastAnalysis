



clear;



%Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('KAS25');

%set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

%turn on or off plotting
analysisParams.showPlots = true;

%Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;
%Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

%set the HRF
load(fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'),'LFContrastAnalysis','subjectHRFs',analysisParams.expSubjID,[analysisParams.expSubjID '_eventGain_results.mat']));
xBase = zeros(1,analysisParams.expLengthTR);
xBase(1:length(results.hrf')) = results.hrf';
analysisParams.HRF.values = xBase;
analysisParams.HRF.timebase =   analysisParams.timebase*1000;
hrfAUC = trapz(analysisParams.HRF.timebase,analysisParams.HRF.values);
analysisParams.HRF.values = analysisParams.HRF.values ./hrfAUC;

[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% angle          = -45;
% minorAxisRatio = 0.19;
% crfAmp         = 5;
% crfExponent    = 1.5;
% crfSemi        = 0.9;
% fullCleanData  = simulateDataFromEllipseParams(analysisParams,angle,minorAxisRatio,'numVoxels',850,...
%                 'crfOffset', 0,'noiseSD',2, 'noiseInverseFrequencyPower', .3,'crfAmp',crfAmp,'crfExponent',crfExponent,'crfSemi',crfSemi);

% Fit IAMP
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, ~, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% Get directon/contrast form of time course and IAMP crf packet pockets.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
directionTimeCoursePacketPocket = {directionTimeCoursePacketPocket{1,:},directionTimeCoursePacketPocket{2,:}};

for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','averageBaseline');
    %[concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end

medianIampParams = iampOBJ.medianParams(concatParams);

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,medianIampParams);


%% Fit error scalar matters
fitErrorScalar  = 1000;


% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket},'fitErrorScalar',fitErrorScalar);

[qcmTCmeanOBJ,qcmTCParams] = fitDirectionModel(analysisParams, 'qcmFit', directionTimeCoursePacketPocket);
qcmTCmeanParams = qcmTCmeanOBJ.medianParams(qcmTCParams);
qcmTCmeanOBJ.paramPrint(qcmTCmeanParams)






