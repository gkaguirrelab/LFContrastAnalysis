function [nrCrfParamsAmpVec, nrCrfParamsExpVec, nrCrfParamsAmpExpVec, qcmCrfMeanParamsVec] = runDirectionModelFits(analysisParams,iampTCPacketPocketBoot,iampParamsbootstrap,iampOBJ)

% Get directon/contrast form of time course and IAMP crf packet pockets.
%
% This conversion is possible because the IAMP packet pocket has meta data
% that we put there to allow exactly this conversion.  That meta data
% encapsulates the key things we need to know about the stimulus obtained
% from the analysis parameters.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTCPacketPocketBoot);

% This puts together pairs of acquistions from the two sessions, so that
% we have one IAMP fit for each pair.  We do this because to fit the
% quadratic model, we need data for all of the color directions together.
%
% NOTE: This bit is very specific to the design of the experiment we are
% currently analyzing, and has to do specifically with the way color
% directions were studied across acquisitions and sessions.
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParamsbootstrap(:,ii),'baselineMethod','makeBaselineZero');
end

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampOBJ.averageParams(concatParams));

%% Fit the direction based models to the mean IAMP beta weights
%
% Fit the CRF with the NR common amplitude -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'talkToMe', false);

% Fit the CRF with the NR common Exponent -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonExp', true, 'talkToMe', false);

% Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonExp', true, 'talkToMe', false);

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket}, 'talkToMe', false);

nrCrfParamsAmpVec = nrCrfOBJ.paramsToVec(nrCrfParamsAmp{1});

nrCrfParamsExpVec = nrCrfOBJ.paramsToVec(nrCrfParamsExp{1});

nrCrfParamsAmpExpVec = nrCrfOBJ.paramsToVec(nrCrfParamsAmpExp{1});

qcmCrfMeanParamsVec = qcmCrfMeanOBJ.paramsToVec(qcmCrfMeanParams{1});



end