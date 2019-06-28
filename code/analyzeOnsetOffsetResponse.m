% Initialize
clear; 

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26_replication');

% SIMULATE MODE
analysisParams.analysisSimulate = false;

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
if analysisParams.analysisSimulate
    analysisParams.numAcquisitions = 10;
    analysisParams.numSessions = 2;
    betaWeights = [repmat(1:-1/5:1/5,1,4), 0]';
    numDirections = 4;
    numContrast = 6; 
    numVoxels = 400;
    [params,fullCleanData] = simulateDataFromExpParams(analysisParams,betaWeights,numDirections,numContrast,numVoxels, 'linDetrending', false);
else
    [fullCleanData, analysisParams] = getTimeCourse(analysisParams);
end

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

%IAMP -- Onset Midoint Offset
[analysisParams, iampOnMidOffPacketPocket, iampOBJ, iampOnMidOffParams, iampOnMidOffResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData, 'modelOnOff', true, 'concatAndFit', true);

%IAMP -- Onset Offset
[~, iampOnOffPacketPocket, ~, iampOnOffParams, iampOnOffResponses, ~] = fit_IAMP(analysisParams,fullCleanData, 'modelOnOff', true, 'concatAndFit', true, 'midpoint', false);

%IAMP -- Block
[~, iampBlockPacketPocket, ~, iampBlockParams, iampBlockResponses, ~] = fit_IAMP(analysisParams,fullCleanData, 'modelOnOff', false, 'concatAndFit', true);


% This puts together pairs of acquistions from the two sessions, so that
% we have one IAMP fit for each pair.  We do this because to fit the
% quadratic model, we need data for all of the color directions together.
%
% NOTE: This bit is very specific to the design of the experiment we are
% currently analyzing, and has to do specifically with the way color
% directions were studied across acquisitions and sessions.

for ii = 1:size(iampOnMidOffParams,2)
    [concatParams{ii},concatBaselineShiftOnMidOff(:,ii)] = iampOBJ.concatenateParams(iampOnMidOffParams(:,ii),'baselineMethod','averageBaseline');
    [concatParams{ii},concatBaselineShiftOnOff(:,ii)] = iampOBJ.concatenateParams(iampOnOffParams(:,ii),'baselineMethod','averageBaseline');
    [concatParams{ii},concatBaselineShiftBlock(:,ii)] = iampOBJ.concatenateParams(iampBlockParams(:,ii),'baselineMethod','averageBaseline');
end

averageIampParams = iampOBJ.averageParams(concatParams);

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,averageIampParams);

% Get the time course prediction from the avarage IAMP params`
iampTCOnMidOff.sessionOne  = iampOBJ.averageParams(iampOnMidOffParams(1,:));
iampTCOnMidOff.sessionTwo  = iampOBJ.averageParams(iampOnMidOffParams(2,:));
iampTCOnMidOff.baseline = concatBaselineShiftOnMidOff;
timeCoursePlot.concatIampOnMidOff = responseFromPacket('meanIAMP', analysisParams, iampTCOnMidOff, iampOnMidOffPacketPocket, 'plotColor', [1 0.0 0.0]);


iampTCOnOff.sessionOne  = iampOBJ.averageParams(iampOnOffParams(1,:));
iampTCOnOff.sessionTwo  = iampOBJ.averageParams(iampOnOffParams(2,:));
iampTCOnOff.baseline = concatBaselineShiftOnOff;
timeCoursePlot.concatIampOnOff = responseFromPacket('meanIAMP', analysisParams, iampTCOnOff, iampOnOffPacketPocket, 'plotColor', [0.0 1.0 0.0]);

iampTCblock.sessionOne  = iampOBJ.averageParams(iampBlockParams(1,:));
iampTCblock.sessionTwo  = iampOBJ.averageParams(iampBlockParams(2,:));
iampTCblock.baseline = concatBaselineShiftBlock;
timeCoursePlot.concatIampBlock = responseFromPacket('meanIAMP', analysisParams, iampTCblock, iampBlockPacketPocket, 'plotColor', [0.0 0.0 1.0]);
% Add clean time
timeCoursePlot.timecourse = rawTC;


%Plot the time course prediction for each run using the different fits to
%the crf

tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, concatBaselineShiftOnMidOff, analysisParams.numSessions);
figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
    [analysisParams.expSubjID,'_OnOffTimeCourse_' analysisParams.sessionNickname '.pdf']);
FigureSave(figNameTc,tcHndl,'pdf');

