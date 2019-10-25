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

%% Get the time series
% from simulation
if analysisParams.analysisSimulate
    analysisParams.numAcquisitions = 10;
    analysisParams.numSessions = 2;
    switch analysisParams.simulationMethod
        case 'IAMP'
            betaWeights = [repmat(1:-1/5:1/5,1,4), 0]';
            numDirections = 4;
            numContrast = 6;
            numVoxels = 400;
            [params,fullCleanData] = simulateDataFromExpParams(analysisParams,betaWeights,numDirections,numContrast,numVoxels, 'linDetrending', false);
        case 'QCM'
            angle = -45;
            minorAxisRatio = 0.19;
            fullCleanData = simulateDataFromEllipseParams(analysisParams,angle,minorAxisRatio,'numVoxels',850,...
                'crfOffset', 0,'noiseSD',2, 'noiseInverseFrequencyPower', .1);
    end
    % From the data
else
    switch analysisParams.preproc
        case 'fmriprep'
            [fullCleanData, analysisParams] = getTimeCourse(analysisParams);
        case 'hcp'
            [fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);
        otherwise
            error('Preprocessing method unknown')
    end
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

% Fit the time course packets with the QCM -- { } is because this expects a cell
directionTimeCoursePacketPocket = {directionTimeCoursePacketPocket{1,:},directionTimeCoursePacketPocket{2,:}};
[qcmOBJ,qcmParams] = fitDirectionModel(analysisParams, 'qcmFit', directionTimeCoursePacketPocket);

% plot the nonlinearity
% if analysisParams.showPlots
%     for ii = 1:length(qcmParams)
%         [nrVals] = plotNakaRushtonFromParams(qcmParams{ii}.crfAmp ,qcmParams{ii}.crfExponent,qcmParams{ii}.crfSemi,...
%             'analysisParams',analysisParams,'plotFunction',true,'savePlot',true);
%     end
% end

% Get the time course predicitions fromt the QCM params fit to the CRF
for ii = 1:length(qcmParams)
    qcmTimeCoursePlot(ii,:) = responseFromPacket('qcmPred', analysisParams, qcmParams{ii}, directionTimeCoursePacketPocket(ii), 'plotColor', [.7, .1, .3]);
end


%% plot the results

figure
counter =0;
for ii = 1:length(qcmParams)
    if ii == 11
        figure
        counter = 0;
    end
    counter = counter +1;
    subplot(3,4,counter)
    hold on
    plot( directionTimeCoursePacketPocket{ii}.response.timebase, directionTimeCoursePacketPocket{ii}.response.values,'k')
    plot(qcmTimeCoursePlot{ii}.timebase,qcmTimeCoursePlot{ii}.values,'color',qcmTimeCoursePlot{ii}.plotColor,'linewidth',2)
    plot(iampResponses{ii}.timebase,iampResponses{ii}.values,'g','linewidth',2)
end

%% create an average QCM from individual run fits
[medianQcmParams] = qcmOBJ.medianParams(qcmParams);
[nrVals] = plotNakaRushtonFromParams(medianQcmParams.crfAmp ,medianQcmParams.crfExponent,medianQcmParams.crfSemi,...
    'analysisParams',analysisParams,'plotFunction',true,'savePlot',true);

% Get the time course predicitions fromt the QCM params fit to the CRF
% Split regressors
medianIampParamsSplit.sessionOne = medianIampParams;
medianIampParamsSplit.sessionOne.paramMainMatrix = [medianIampParams.paramMainMatrix(1:20);medianIampParams.paramMainMatrix(end)];
medianIampParamsSplit.sessionOne.matrixRows      = length(medianIampParamsSplit.sessionOne.paramMainMatrix);

medianIampParamsSplit.sessionTwo = medianIampParams;
medianIampParamsSplit.sessionTwo.paramMainMatrix = [medianIampParams.paramMainMatrix(21:40);medianIampParams.paramMainMatrix(end)];
medianIampParamsSplit.sessionTwo.matrixRows      = length(medianIampParamsSplit.sessionTwo.paramMainMatrix);

for ii= 1:length(directionTimeCoursePacketPocket)
    qcmMedianTimeCoursePlot(ii,:)   = responseFromPacket('qcmPred', analysisParams, medianQcmParams, directionTimeCoursePacketPocket(ii), 'plotColor', [.3, .1, .7]);
end
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
iampMedianTimeCoursePlot  = responseFromPacket('meanIAMP', analysisParams, medianIampParamsSplit, directionTimeCoursePacketPocket, 'plotColor', [.4, .5, .35]);
iampMedianTimeCoursePlot  = {iampMedianTimeCoursePlot{1,:},iampMedianTimeCoursePlot{2,:}};


%% plot the results
figure
counter =0;
for ii = 1:length(qcmParams)
    if ii == 11
        figure
        counter = 0;
    end
    counter = counter +1;
    subplot(3,4,counter)
    hold on
    plot( directionTimeCoursePacketPocket{ii}.response.timebase, directionTimeCoursePacketPocket{ii}.response.values,'k')
    plot(qcmMedianTimeCoursePlot{ii}.timebase,qcmMedianTimeCoursePlot{ii}.values,'color',qcmMedianTimeCoursePlot{ii}.plotColor,'linewidth',2)
    plot(iampMedianTimeCoursePlot{ii}.timebase,iampMedianTimeCoursePlot{ii}.values,'color',iampMedianTimeCoursePlot{ii}.plotColor,'linewidth',2)
end