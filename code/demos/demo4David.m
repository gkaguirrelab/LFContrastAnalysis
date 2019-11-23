%% This is a demo to help us fit the QCM directly to the time course
%

%% Initialize
clear; close all

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26');

% Set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

% Turn on or off plotting
analysisParams.showPlots = true;

% Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;

% Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

% Using the canonical HRF until I meet with Geoff for fitted HRFs
analysisParams.HRF = generateHRFKernel(6,12,10,analysisParams.timebase*1000);

% Get the data
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

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
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData,'highpass',false);

% reshape
iampResponses = {iampResponses{1,:},iampResponses{2,:}};

% Get directon/contrast form of time course and IAMP crf packet pockets.
%
% This conversion is possible because the IAMP packet pocket has meta data
% that we put there to allow exactly this conversion.  That meta data
% encapsulates the key things we need to know about the stimulus obtained
% from the analysis parameters.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);

% ###### FIX ###################
% remove subraction of the baseline
% ##############################
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','averageBaseline');
    %[concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end

%% Get median IAMP parameters
medianIampParams = iampOBJ.medianParams(concatParams);
directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,medianIampParams);

%% Which run to look at in detail
runIdx = 3;
directionTimeCoursePacket = directionTimeCoursePacketPocket{runIdx};

%% Fit error scalar matters
fitErrorScalar  = 10000;

%% Fit QCM model to the IAMP parameters 
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit',{directionCrfMeanPacket},'fitErrorScalar',fitErrorScalar);
qcmTimeCourseCrf = responseFromPacket('qcmPred', analysisParams, qcmCrfMeanParams{1},{directionTimeCoursePacket});
fQcmTimeCourseCrf = qcmCrfMeanOBJ.fitError(qcmCrfMeanOBJ.paramsToVec(qcmCrfMeanParams{1}),directionTimeCoursePacket,'fitErrorScalar',fitErrorScalar);

%% Fit the time course packets with the QCM -- { } is because this expects a cell
directionTimeCoursePacketPocket = {directionTimeCoursePacketPocket{1,:},directionTimeCoursePacketPocket{2,:}};

% Unseeded fit and time course predictions
[qcmOBJ,qcmParamsUnseeded] = fitDirectionModel(analysisParams, 'qcmFit',{directionTimeCoursePacket},'fitErrorScalar',fitErrorScalar);
qcmTimeCourseUnseeded = responseFromPacket('qcmPred', analysisParams, qcmParamsUnseeded{1},{directionTimeCoursePacket});
fQcmTimeCourseUnseeded = qcmOBJ.fitError(qcmOBJ.paramsToVec(qcmParamsUnseeded{1}),directionTimeCoursePacket,'fitErrorScalar',fitErrorScalar);

% Seeded fit and time course predictions
[qcmOBJ,qcmParamsSeeded] = fitDirectionModel(analysisParams, 'qcmFit',{directionTimeCoursePacket},'initialParams',qcmCrfMeanParams,'fitErrorScalar',fitErrorScalar);
qcmTimeCourseSeeded = responseFromPacket('qcmPred', analysisParams, qcmParamsSeeded{1},{directionTimeCoursePacket});
fQcmTimeCourseSeeded = qcmOBJ.fitError(qcmOBJ.paramsToVec(qcmParamsSeeded{1}),directionTimeCoursePacket,'fitErrorScalar',fitErrorScalar);

% Plot the fits 
figure;
hold on
plot(qcmTimeCourseUnseeded{1}.timebase,qcmTimeCourseCrf{1}.values,'g','LineWidth',4);
plot(qcmTimeCourseUnseeded{1}.timebase,qcmTimeCourseUnseeded{1}.values,'r','LineWidth',4);
plot(qcmTimeCourseSeeded{1}.timebase,qcmTimeCourseSeeded{1}.values,'b','LineWidth',2);
plot(directionTimeCoursePacket.response.timebase,directionTimeCoursePacket.response.values,'k','LineWidth',2);
legend('CRF Fit','Unseeded Fit','Seeded Fit','Time Course')
xlabel('Time (s)')

% Report error
title(sprintf('Fit error CRF: %g; Unseeded: %g; Seeded: %g',fQcmTimeCourseCrf,fQcmTimeCourseUnseeded,fQcmTimeCourseSeeded));





