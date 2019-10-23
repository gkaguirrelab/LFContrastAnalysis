function [qcmParams,meanRsquaredIAMP,stdRsquaredIAMP, meanRsquaredQCM,stdRsquaredQCM] = fitQCMtoVoxel(analysisParams,voxelTimeSeries)
clear defaultParamsInfo
defaultParamsInfo.noOffset = false;
fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);

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
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,voxelTimeSeries);



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
[qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket},'talkToMe',false);

for ii = 1: size(directionTimeCoursePacketPocket,1)
    if ii == 1
        averageIampParamsSplit = averageIampParams;
        averageIampParamsSplit.paramMainMatrix = [averageIampParams.paramMainMatrix(1:20); averageIampParams.paramMainMatrix(end)];
        averageIampParamsSplit.matrixRows = 21;
    elseif ii == 2
        averageIampParamsSplit = averageIampParams;
        averageIampParamsSplit.paramMainMatrix = [averageIampParams.paramMainMatrix(21:40); averageIampParams.paramMainMatrix(end)];
        averageIampParamsSplit.matrixRows = 21;
    else
        error('not coded up for more than 2 sessions')
    end
    for jj = 1: size(directionTimeCoursePacketPocket,2)
        
        % compute IAMP preditciotn time course and R^2
        
        iampPred = iampOBJ.computeResponse(averageIampParamsSplit,iampTimeCoursePacketPocket{ii,jj}.stimulus,...
            iampTimeCoursePacketPocket{ii,jj}.kernel);
        corrVecIAMP = [iampPred.values',iampTimeCoursePacketPocket{ii,jj}.response.values'];
        corrValsIAMP = corr(corrVecIAMP);
        rSquaredIAMPAllRuns(ii,jj) = corrValsIAMP(1,2).^2;
        
        % compute QCM preditciotn time course and R^2
        qcmPred = fitOBJ.computeResponse(qcmCrfMeanParams{1},directionTimeCoursePacketPocket{ii,jj}.stimulus,...
            directionTimeCoursePacketPocket{ii,jj}.kernel);
        corrVecQCM = [qcmPred.values',directionTimeCoursePacketPocket{ii,jj}.response.values'];
        corrValsQCM = corr(corrVecQCM);
        rSquaredQCMAllRuns(ii,jj) = corrValsQCM(1,2).^2;
        
    end
end
meanRsquaredIAMP = mean(rSquaredIAMPAllRuns(:));
stdRsquaredIAMP  = std(rSquaredIAMPAllRuns(:));
meanRsquaredQCM = mean(rSquaredQCMAllRuns(:));
stdRsquaredQCM  = std(rSquaredQCMAllRuns(:));
qcmParams = qcmCrfMeanOBJ.paramsToVec(qcmCrfMeanParams{1});
