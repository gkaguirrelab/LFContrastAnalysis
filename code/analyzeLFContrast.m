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
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses] = fit_IAMP(analysisParams,fullCleanData);

% Get directon/contrast form of time course and IAMP crf packet pockets
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);

% Seperate out the fits per session to take the average 
for ii = 1:analysisParams.numAcquisitions
    concatParams{ii} = iampOBJ.concatenateParams(iampParams(:,ii));
end

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampOBJ.averageParams(concatParams));

%% Fit the direction based models to the mean IAMP beta weights 

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});

% Fit the CRF with the NR common amplitude -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true);

% Fit the CRF with the NR common amplitude and semisaturation  -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmpSemi] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonSemi', true);

% Fit the CRF with the NR common amplitude, semisaturation, and exponent  -- { } is because this expects a cell
[nrCrfOBJ,nrCrfParamsAmpSemiExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonSemi', true, 'commonExp', true);


[responses] =  responseFromPacket(qcmCrfMeanOBJ, analysisParams, qcmCrfMeanParams, {directionCrfMeanPacket}, 'upsampleCrf', true);
% Upsample the NR repsonses 







%% Compute responses to the direction stimuli from each run with the fit from above


[responses, plotPocket] =  responseFromPacket(obj, params, packetPocket)

 
%% Plot the CRF from the IAMP, QCM, and  fits
% nrParams = plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit);
% 

% %Plot the time course prediction for each run using the different fits to
% %the crf


% % Plot isoresponce contour
% thresholds = [0.10, 0.2, 0.3];
% colors     = [0.5,0.0,0.0; 0.5,0.5,0.0; 0.0,0.5,0.5;];
% [hdl] = plotIsoresponse(analysisParams,meanIAMPBetas,paramsQCMFit,thresholds,nrParams,colors);
% 
% % Use QCM fit to IAMP to predict timecourse.
% plotQCMtimecourse(paramsFitIAMP,packetPocket,meanIAMPBetas,analysisParams,fitResponseStructQCM,paramsQCMFit.crfOffset);
