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
% Load the subject relevant info
analysisParams = getSubjectParams('AP26');

analysisParams.preproc = 'hcp';

theDimension = 2;

%turn on or off plotting
analysisParams.showPlots = true;
qcmColor  = [0.4078, 0.2784, 0.5765];
lcmColor  = [90,171,97]./255;

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
LCMObj = tfeLCMDirection('verbosity','none','dimension',theDimension);

% Fit the LCM packet
[fitLCMParams,fVal,fitResponseStructLCM] = LCMObj.fitResponse(thePacketLcm);
fprintf('\nLCM parameters from fit:\n');
LCMObj.paramPrint(fitLCMParams)
% Fit the time course with the QCM -- { } is because this expects a cell
[qcmTcOBJ,qcmParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000,'talkToMe',false);


%% Get the time course response
LCMResponseStruct = LCMObj.computeResponse(fitLCMParams,thePacketLcm.stimulus,thePacketLcm.kernel,'AddNoise',false);
LCMResponseStruct.plotColor = lcmColor;
[timeCoursePlot.lcm]   = chopUpTimeCourse(LCMResponseStruct,20);

% Get the time course predicitions fromt the QCM params fit to the CRF
qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmParams{1}, timeCoursePacket, 'plotColor', qcmColor);
[timeCoursePlot.qcm] = chopUpTimeCourse(qcmTimeCourse{1},20);

% Add clean time
theTimeCourse = {theFullPacket.response};
theTimeCourse{1}.plotColor =[0, 0, 0];
[timeCoursePlot.timecourse] = chopUpTimeCourse(theTimeCourse{1},20);


%% Plot it
figSavePath =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID);

tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, zeros(20,1), 20);
% Save it
% Plot configuration
set(tcHndl, 'Renderer', 'Painters');
figureSizeInches = [20 15];
set(tcHndl, 'PaperUnits', 'inches');
set(tcHndl, 'PaperSize',figureSizeInches);
set(tcHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
% Full file name
figNameTc =  fullfile(figSavePath,[analysisParams.expSubjID,'_LCM_TimeCourse_' analysisParams.sessionNickname...
    '_' analysisParams.preproc '.pdf']);
% Save it
print(tcHndl, figNameTc, '-dpdf', '-r300');
%% Get and plot isoresponse contour
ellipseHndl = figure; hold on
[isoContrastFit,unitContrastResponseFit,angleSupport] = LCMObj.getIsoContrast(fitLCMParams);
plot(isoContrastFit.*cosd(angleSupport),isoContrastFit.*sind(angleSupport),'r','LineWidth',2, 'color', lcmColor);
axis('square');
xlim([-2 2]); ylim([-2 2]);
xlabel('Cone 1 Contrast');
ylabel('Cone 2 Contrast');
title('IsoContrast');
%% Ellipse Figure
nQCMPoints = 200;
% Calculate the Minv matrix to tranform a unit circle to the ellipse and do it
[~,Minv,~] = EllipsoidMatricesGenerate([1 1./qcmParams{1}.Qvec(1) qcmParams{1}.Qvec(2)]','dimension',2);
circlePoints = UnitCircleGenerate(nQCMPoints);
ellipsePoints = Minv*circlePoints;
% get current axes
axh = gca;

% plot ellipse
line(ellipsePoints(1,:),ellipsePoints(2,:),'color', qcmColor, 'LineWidth', 2);

% plot axes
line([-2 2], [0 0], 'Color', [.3 .3 .3], 'LineStyle', ':','LineWidth', 2);
line([0 0], [-2 2], 'Color', [.3 .3 .3], 'LineStyle', ':','LineWidth', 2);

legend('LCM','QCM')

set(ellipseHndl, 'Renderer', 'Painters');
figureSizeInches = [7 6];
set(ellipseHndl, 'PaperUnits', 'inches');
set(ellipseHndl, 'PaperSize',figureSizeInches);
set(ellipseHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
figNameEllipseNonlin = fullfile(figSavePath,[analysisParams.expSubjID,'_LCM_Ellipse_' ...
    analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
print(ellipseHndl, figNameEllipseNonlin, '-dpdf', '-r300');


