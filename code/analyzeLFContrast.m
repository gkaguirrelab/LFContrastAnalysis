% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('LZ23');

% Clip fisrt 2 TRs from time series?
% if no clipping then put 0;
analysisParams.numClipFramesStart = 0;
analysisParams.numClipFramesEnd   = 2;

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

% Fit IAMP 
% 
% Fit IAMP to each constructed packet and create packetPocket cell array of
% all the fit packets.
%     packetPocket - Meta data of packePocket contains the direction/contrast form of the same packet.
%     iampOBJ - the tfe IAMP object
%     iampParams - cell array of iampParams for each object   
[analysisParams,iampTimeCoursePacketPocket,iampOBJ,iampParams] = fit_IAMP(analysisParams,fullCleanData);

% Get directon/contrast form of time course and IAMP crf packet pockets
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);

% Seperate out the fits per session to take the average 
for ii = 1:analysisParams.numAcquisitions
    concatParams{ii} = iampOBJ.concatenateParams(iampParams(:,ii));
end

directionCrfMeanPacketPocket = makeDirectionCrfPacketPocket(analysisParams,iampOBJ.averageParams(concatParams));
% % Fit the direction based models
% % 
% % Here is an example for QCM
% [qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel('qcmFit',analysisParams,directionCrfMeanPacketPocket);

% 
% % Plot the CRF from the IAMP and QCM fits
% nrParams = plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit);
% 
% % Plot isoresponce contour
% thresholds = [0.10, 0.2, 0.3];
% colors     = [0.5,0.0,0.0; 0.5,0.5,0.0; 0.0,0.5,0.5;];
% [hdl] = plotIsoresponse(analysisParams,meanIAMPBetas,paramsQCMFit,thresholds,nrParams,colors);
% 
% % Use QCM fit to IAMP to predict timecourse.
% plotQCMtimecourse(paramsFitIAMP,packetPocket,meanIAMPBetas,analysisParams,fitResponseStructQCM,paramsQCMFit.crfOffset);
