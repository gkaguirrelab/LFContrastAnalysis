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
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% Get directon/contrast form of time course and IAMP crf packet pockets.
%
% This conversion is possible because the IAMP packet pocket has meta data
% that we put there to allow exactly this conversion.  That meta data
% encapsulates the key things we need to know about the stimulus obtained
% from the analysis parameters.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);

% This puts together pairs of acquistions from the two sessions, so that
% we have one IAMP fit for each pair.  We do this because to fit the
% quadratic model, we need data for all of the color directions together.
%
% NOTE: This bit is very specific to the design of the experiment we are
% currently analyzing, and has to do specifically with the way color
% directions were studied across acquisitions and sessions.


% ###### FIX ###################
% remove subraction of the baseline
% ##############################
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','averageBaseline');
    %[concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end

averageIampParams = iampOBJ.averageParams(concatParams);

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,averageIampParams);

%% Fit the direction based models to the mean IAMP beta weights
%



% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});


[qcmTimeCourseOBJ,qcmTimeCourseParams] = fitDirectionModel(analysisParams, 'qcmFit', directionTimeCoursePacketPocket);


avgQCMparams = qcmTimeCourseOBJ.averageParams(qcmTimeCourseParams);

% Get the time course predicitions of the CRF params

% Get the time course predicitions fromt the QCM params fit to the CRF
timeCoursePlot.QCM = responseFromPacket('qcmPred', analysisParams, qcmCrfMeanParams{1}, directionTimeCoursePacketPocket, 'plotColor', [0, 1, 0]);
timeCoursePlot.averageQCM = responseFromPacket('qcmPred', analysisParams, avgQCMparams, directionTimeCoursePacketPocket, 'plotColor', [1, 0, 0]);


% Get the predictions from individual IAMP params 
for ii = 1:size(iampParams,1)
    for jj = 1:size(iampParams,2)
        timeCoursePlot.iamp{ii,jj} = responseFromPacket('IAMP', analysisParams, iampParams{ii,jj}, iampTimeCoursePacketPocket{ii,jj}, 'plotColor', [0.5 0.2 0]);
    end
end

%% FIX THIS
% Get the time course prediction from the avarage IAMP params
iampParamsTC.sessionOne  = iampOBJ.averageParams(iampParams(1,:));
iampParamsTC.sessionTwo  = iampOBJ.averageParams(iampParams(2,:));
iampParamsTC.baseline = concatBaselineShift;
timeCoursePlot.meanIamp = responseFromPacket('meanIAMP', analysisParams, iampParamsTC, iampTimeCoursePacketPocket, 'plotColor', [0.0 0.35 0.9]);

% Add clean time
timeCoursePlot.timecourse = rawTC;


%Plot the time course prediction for each run using the different fits to
%the crf

tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, concatBaselineShift, analysisParams.numSessions*analysisParams.numAcquisitions);
figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
    [analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname '.pdf']);
FigureSave(figNameTc,tcHndl,'pdf');

