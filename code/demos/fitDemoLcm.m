%% DEMO for fitting the LCM with data from LFContrast
%
% Description:
%   This script synthesizes data for the LCM model and then fits it,
%   to ensure that we can get back what we put in.  Etc.
%
%   The parameters used to synthesize the data were just pulled out of a
%   hat, so there is no particular expected shape the resulting isoresponse
%   contour.

% History:
%   03/01/21  dhb       Wrote from QCM test version.
%   03/07/21  mab       Integrated LCM with real data

%% Initialize
clear all; rng(0);
% Fit error scalar.  Big value seems
% to work better here.
fitErrorScalar = 10000;

% Load the subject relevant info
analysisParams = getSubjectParams('KAS25_replication');

analysisParams.preproc = 'hcp';

theDimension = 2;

%turn on or off plotting
analysisParams.showPlots = true;
qcmColor  = [0.4078, 0.2784, 0.5765];
lcmColorBH  = [90,171,97]./255;
lcmColorKim = [171, 109, 89]./255;

%% Set up the LCM

summationExponent = 1;
startCenter = 0;
% This implements our version of the Brouwer and Heeger model.
% Six channels with cos^2 sensitivity, linear combination.
%
% Channel properties
nChannelsBH = 6;
channelExponentBH = 2;
LCMObjBH = tfeLCMDirection('verbosity','none','dimension',theDimension, ...
    'nChannels',nChannelsBH,'channelExponent',channelExponentBH,'summationExponent',...
    summationExponent,'startCenter',startCenter);
% This implements the Kim et al. variant of
% the Brouwer and Heeger model.
% Eight channels with cos^6 sensitivity, linear combination.
%
% Channel properties
nChannelsKim = 8;
channelExponentKim = 6;
LCMObjKim = tfeLCMDirection('verbosity','none','dimension',theDimension, ...
    'nChannels',nChannelsKim,'channelExponent',channelExponentKim,'summationExponent',...
    summationExponent,'startCenter',startCenter);

%% Load the relevant data (SDM, HRF, TC)

%set the HRF
[analysisParams] = loadHRF(analysisParams);

% Load the time course
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Get a packet for each run (1-20)
[analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData);

% Concatenate Packet
[analysisParams, theFullPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);

% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCoursePacket = makeDirectionTimeCoursePacketPocket({theFullPacket});

thePacketLcm = timeCoursePacket{1};

%% FIT THE TIME COURSE
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');
% Get the LCM fut object


% Fit the LCM packet BH
[fitLCMParamsBH,~,~] =LCMObjBH.fitResponse(thePacketLcm,'fitErrorScalar',fitErrorScalar,...
    'maxIter',3000,'maxFunEval',3000);
fprintf('\nLCM parameters from fit:\n');
LCMObjBH.paramPrint(fitLCMParamsBH)

% Fit the LCM packet Kim
[fitLCMParamsKim,~,~] =LCMObjKim.fitResponse(thePacketLcm,'fitErrorScalar',fitErrorScalar,...
    'maxIter',3000,'maxFunEval',3000);
fprintf('\nLCM parameters from fit:\n');
LCMObjKim.paramPrint(fitLCMParamsKim)

LCMObjBH.criterionResp = 0.3;
isoContrastLcmBH = LCMObjBH.getIsoContrast(fitLCMParamsBH);

LCMObjKim.criterionResp = 0.3;
isoContrastLcmKim = LCMObjKim.getIsoContrast(fitLCMParamsKim);

% Fit the time course with the QCM -- { } is because this expects a cell
[qcmTcOBJ,qcmParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000,'talkToMe',false);

%% Generate an ellipsoidal isoresponse contour
criterionRespQCM =1;
ellipticalIsoContrast = tfeEllipticalIsoContrast(qcmParams{1}.Qvec(2),qcmParams{1}.Qvec(1),LCMObjBH.angleSupport,criterionRespQCM);

% Find scalar
scaleFactorBH = isoContrastLcmBH'\ellipticalIsoContrast';
scaleFactorKim = isoContrastLcmKim'\ellipticalIsoContrast';

scaledIsoContrastLcmBH = isoContrastLcmBH .* scaleFactorBH;
scaledIsoContrastLcmKim = isoContrastLcmKim .* scaleFactorKim;

% %% Get the time course response
% LCMResponseStruct = LCMObj.computeResponse(fitLCMParams,thePacketLcm.stimulus,thePacketLcm.kernel,'AddNoise',false);
% LCMResponseStruct.plotColor = lcmColor;
% [timeCoursePlot.lcm]   = chopUpTimeCourse(LCMResponseStruct,20);
% 
% % Get the time course predicitions fromt the QCM params fit to the CRF
% qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmParams{1}, timeCoursePacket, 'plotColor', qcmColor);
% [timeCoursePlot.qcm] = chopUpTimeCourse(qcmTimeCourse{1},20);
% 
% % Add clean time
% theTimeCourse = {theFullPacket.response};
% theTimeCourse{1}.plotColor =[0, 0, 0];
% [timeCoursePlot.timecourse] = chopUpTimeCourse(theTimeCourse{1},20);

% 
% %% Plot it time course
% figSavePath =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID);
% 
% tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, zeros(20,1), 20);
% % Save it
% % Plot configuration
% set(tcHndl, 'Renderer', 'Painters');
% figureSizeInches = [20 15];
% set(tcHndl, 'PaperUnits', 'inches');
% set(tcHndl, 'PaperSize',figureSizeInches);
% set(tcHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
% % Full file name
% figNameTc =  fullfile(figSavePath,[analysisParams.expSubjID,'_LCM_TimeCourse_' analysisParams.sessionNickname...
%     '_' analysisParams.preproc '.pdf']);
% % Save it
% print(tcHndl, figNameTc, '-dpdf', '-r300');


%% Get and plot isoresponse contour
ellipseHndl = figure; hold on

%plot axes
line([-2 2], [0 0], 'Color', [.3 .3 .3], 'LineStyle', ':','LineWidth', 2);
line([0 0], [-2 2], 'Color', [.3 .3 .3], 'LineStyle', ':','LineWidth', 2);

y1 = plot(ellipticalIsoContrast.*cosd(LCMObjBH.angleSupport),ellipticalIsoContrast.*sind(LCMObjBH.angleSupport),'Color',qcmColor,'LineWidth',2);
y2 = plot(scaledIsoContrastLcmBH.*cosd(LCMObjBH.angleSupport),scaledIsoContrastLcmBH.*sind(LCMObjBH.angleSupport),'Color',lcmColorBH,'LineWidth',2);
y3 = plot(scaledIsoContrastLcmKim.*cosd(LCMObjKim.angleSupport),scaledIsoContrastLcmKim.*sind(LCMObjKim.angleSupport),'Color',lcmColorKim,'LineWidth',2);

axis('square');
xlim([-1 1]); ylim([-1 1]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');

legend([y1,y2,y3],{'QCM', 'LCM BH', 'LCM Kim'});
set(ellipseHndl, 'Renderer', 'Painters');
figureSizeInches = [7 6];
set(ellipseHndl, 'PaperUnits', 'inches');
set(ellipseHndl, 'PaperSize',figureSizeInches);
set(ellipseHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
figNameEllipseNonlin = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
    [analysisParams.expSubjID,'_LCM_Ellipse_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
print(ellipseHndl, figNameEllipseNonlin, '-dpdf', '-r300');


