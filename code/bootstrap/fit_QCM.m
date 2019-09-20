function [qcmParams, nrVals] = fit_QCM(analysisParams,iampParams)

iampOBJ = tfeIAMP('verbosity','none');

for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','averageBaseline');
    %[concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end
averageIampParams = iampOBJ.averageParams(concatParams);

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,averageIampParams);

%% Fit the direction based models to the mean IAMP beta weights
%

% Fit the CRF -- { } is because this expects a cell
%[nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket});

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket},'talkToMe',false);
qcmParams = qcmCrfMeanOBJ.paramsToVec(qcmCrfMeanParams{1});

% Do some plotting of these fits
nrVals = plotNakaRushtonFromParams(qcmCrfMeanParams{1}.crfAmp ,qcmCrfMeanParams{1}.crfExponent,qcmCrfMeanParams{1}.crfSemi,...
    'analysisParams',analysisParams,'plotFunction',false,'savePlot',false);

end
