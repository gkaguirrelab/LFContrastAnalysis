%% This is a demo to help us fit the QCM directly to the time course
%

%% Initialize
clear;

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26');

% set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

% turn on or off plotting
analysisParams.showPlots = true;

% Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;

% Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

% using the canonical HRF until I meet with Geoff for fitted HRFs
analysisParams.HRF = generateHRFKernel(6,12,10,analysisParams.timebase*1000);

% get the data
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

medianIampParams = iampOBJ.medianParams(concatParams);
directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,medianIampParams);

[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit',{ directionCrfMeanPacket });

% Fit the time course packets with the QCM -- { } is because this expects a cell
directionTimeCoursePacketPocket = {directionTimeCoursePacketPocket{1,:},directionTimeCoursePacketPocket{2,:}};


%% Fitting with QCM for example run that show a difference in fit when seeded
runIdx = 3;
% unseeded fit 
[qcmOBJ,qcmParamsUnseeded] = fitDirectionModel(analysisParams, 'qcmFit', directionTimeCoursePacketPocket(runIdx));

% seeded fit
[qcmOBJ,qcmParamsSeeded] = fitDirectionModel(analysisParams, 'qcmFit', directionTimeCoursePacketPocket(runIdx),'initialParams',qcmCrfMeanParams);

% generate timecourse predictions unseeded
qcmTimeCourseUnseeded = responseFromPacket('qcmPred', analysisParams, qcmParamsUnseeded{1}, directionTimeCoursePacketPocket(runIdx));
% generate timecourse predictions seeded
qcmTimeCourseSeeded = responseFromPacket('qcmPred', analysisParams, qcmParamsSeeded{1}, directionTimeCoursePacketPocket(runIdx));

% plot the fits 
figure;
hold on 
plot(qcmTimeCourseUnseeded{1}.timebase,qcmTimeCourseUnseeded{1}.values,'r')
plot(qcmTimeCourseSeeded{1}.timebase,qcmTimeCourseSeeded{1}.values,'b')
plot(directionTimeCoursePacketPocket{runIdx}.response.timebase,directionTimeCoursePacketPocket{runIdx}.response.values,'k')
legend('Unseeded Fit','Seeded Fit','Time Course')
xlabel('Time (s)')





