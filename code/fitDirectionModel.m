
%% Fit the tfeQCM to the IAMP beta weights
% allow QCM to fit the offset
clear defaultParamsInfo
defaultParamsInfo.noOffset = false;
[paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket,'defaultParamsInfo',defaultParamsInfo);
fprintf('\nModel parameter from fits from tfeQCM to IAMP betas:\n');
temporalFitQCM.paramPrint(paramsQCMFit)

%% Fit the tfeQCMDirections to the IAMP beta weights
% create the packet
stim = kron(analysisParams.directionCoding(1:analysisParams.theDimension,:).*analysisParams.maxContrastPerDir,analysisParams.contrastCoding);
stim = [stim, [0;0]];
[qcmDirStimDirections,qcmDirStimContrasts] = tfeQCMStimuliToDirectionsContrasts(stim,'precision',4);
qcmDirPacket.stimulus.values   = [qcmDirStimDirections; qcmDirStimContrasts];
qcmDirPacket.stimulus.timebase = 1:size(qcmDirPacket.stimulus.values,2);
qcmDirPacket.response.values   = meanIAMPBetas';
qcmDirPacket.response.timebase = 1:length(meanIAMPBetas);
qcmDirPacket.kernel            = [];
qcmDirPacket.metaData          = [];

% Create the tfeQCMDirection object
clear defaultParamsInfo
defaultParamsInfo.noOffset = false;
QCMDirectionObj = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);

% Fit the packet
[fitParams.QCMDirParams,fVal,fitQCMDirectionResponseStruct] = QCMDirectionObj.fitResponse(qcmDirPacket,'defaultParamsInfo',defaultParamsInfo);
fprintf('\nQCMDirection parameters from direction fit to IAMP betas:\n');
QCMDirectionObj.paramPrint(fitParams.QCMDirParams)

% Create the packet for NR fits
nrDirPacket = qcmDirPacket;
[uniqueDirections,directionIndices] = tfeQCMParseDirections(qcmDirStimDirections,'precision',4);

%% Fit the NRDirections to the the IAMP beta weights WITH COMMON OFFSET
% Create the tfeNakaRushtonDirection object
commonOffsetNRObj = tfeNakaRushtonDirection(uniqueDirections, ...
    'lockOffsetToZero',false,'commonAmp',false,'commonSemi',false,'commonExp',false,'commonOffset',true);

% Fit the packet
[fitParams.NRDirParamsOff,~,objFitResponses] = commonOffsetNRObj.fitResponse(nrDirPacket);
fprintf('\nNRDirection parameters from fit to IAMP betas with common offset:\n');
commonOffsetNRObj.paramPrint(fitParams.NRDirParamsOff)

%% Fit the NRDirections to the the IAMP beta weights WITH COMMON OFFSET & AMP
% Create the tfeNakaRushtonDirection object
commonAmpNRObj = tfeNakaRushtonDirection(uniqueDirections, ...
    'lockOffsetToZero',false,'commonAmp',true,'commonSemi',false,'commonExp',false,'commonOffset',true);

% Fit the packet
[fitParams.NRDirParamsOffAmp,~,objFitResponses] = commonAmpNRObj.fitResponse(nrDirPacket);
fprintf('\nNRDirection parameters from fit to IAMP betas with common offset & amplitude:\n');
commonAmpNRObj.paramPrint(fitParams.NRDirParamsOffAmp)

%% Fit the NRDirections to the the IAMP beta weights WITH COMMON OFFSET, AMP, & EXP
% Create the tfeNakaRushtonDirection object
commonAmpExpNRObj = tfeNakaRushtonDirection(uniqueDirections, ...
    'lockOffsetToZero',false,'commonAmp',true,'commonSemi',false,'commonExp',true,'commonOffset',true);

% Fit the packet
[fitParams.NRDirParamsOffAmpExp,~,objFitResponses] = commonAmpExpNRObj.fitResponse(nrDirPacket);
fprintf('\nNRDirection parameters from fit to IAMP betas with common offset & amplitude & exponent:\n');
commonAmpExpNRObj.paramPrint(fitParams.NRDirParamsOffAmpExp)