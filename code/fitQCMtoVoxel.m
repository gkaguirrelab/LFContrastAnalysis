function [qcmParams,meanRsquaredIAMP,stdRsquaredIAMP, meanRsquaredQCM,stdRsquaredQCM] = fitQCMtoVoxel(fitOBJ, analysisParams,voxelTimeSeries)


%set the HRF
[analysisParams] = loadHRF(analysisParams);

if analysisParams.highpass
    analysisParams.HRF.values = highpass(analysisParams.HRF.values ,5/288,1/.8);
end

% Load the time course
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Get a packet for each run (1-20)
[analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData,'highpass',analysisParams.highpass);

%% Generate a cell array of concat train packets and a corresponding cell
%  array of concat test packets
[theTestPackets, theTrainPackets,leaveOutPairs] = concatPackets_crossVal(analysisParams, iampTimeCoursePacketPocket);

% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCourseTrainPackets = makeDirectionTimeCoursePacketPocket(theTrainPackets);
timeCourseTestPackets = makeDirectionTimeCoursePacketPocket(theTestPackets);
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');


timeCoursePlot.qcm = [];
timeCoursePlot.IAMP = [];
timeCoursePlot.timecourse = [];
rSquaredQcm = [];
rSquaredIamp = [];
count = 0;
for ii = 1:length(theTrainPackets)
    
    %% FIT THE TIME COURSE
    trainPacket = theTrainPackets{ii};
    timeCoursePacket = timeCourseTrainPackets{ii};
    tcTestPacket = timeCourseTestPackets{ii};
    iampTestPacket = theTestPackets{ii};
    
    % fit the IAMP model
    defaultParamsInfo.nInstances = size(trainPacket.stimulus.values,1);
    [iampParams,fVal,iampResponses] = iampOBJ.fitResponse(trainPacket,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    
    
    % Fit the time course with the QCM -- { } is because this expects a cell
    [qcmTcOBJ,qcmTcParams] = fitDirectionModel(analysisParams, 'qcmFit', {timeCoursePacket},'fitErrorScalar',1000,'talkToMe',false);
    
    % Get the time course predicitions fromt the QCM params fit to the CRF
    qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmTcParams{1}, {tcTestPacket}, 'plotColor', qcmColor);
    qcmChopped = chopUpTimeCourse(qcmTimeCourse{1},2);
    [timeCoursePlot.qcm] = [timeCoursePlot.qcm,qcmChopped];
    
    iampTimeCourse = responseFromPacket('IAMP', analysisParams, iampParams, iampTestPacket, 'plotColor', iampColor);
    iampChopped = chopUpTimeCourse(iampTimeCourse,2);
    [timeCoursePlot.IAMP] = [timeCoursePlot.IAMP, iampChopped];
    
    theTimeCourse= tcTestPacket.response;
    theTimeCourse.plotColor =[0, 0, 0];
    timeCourseChopped = chopUpTimeCourse(theTimeCourse,2);
    [timeCoursePlot.timecourse] = [timeCoursePlot.timecourse, timeCourseChopped];
    
    %% Calc R squared
    
    for jj = 1:length(timeCourseChopped)
        count = count +1;
        qcmCorrVec =  [timeCourseChopped{jj}.values',qcmChopped{jj}.values'];
        iampCorrVec =  [timeCourseChopped{jj}.values',iampChopped{jj}.values'];
        
        qcmCorrVals = corrcoef(qcmCorrVec(:,1),qcmCorrVec(:,2),'rows','complete').^2;
        iampCorrVals = corrcoef(iampCorrVec(:,1),iampCorrVec(:,2),'rows','complete').^2;
        
        rSquaredQcm(count) = qcmCorrVals(1,2);
        rSquaredIamp(count) = iampCorrVals(1,2);
    end
    
end


meanRsquaredIAMP = mean(rSquaredIAMPAllRuns(:));
stdRsquaredIAMP  = std(rSquaredIAMPAllRuns(:));
meanRsquaredQCM = mean(rSquaredQCMAllRuns(:));
stdRsquaredQCM  = std(rSquaredQCMAllRuns(:));
qcmParams = qcmCrfMeanOBJ.paramsToVec(qcmCrfMeanParams{1});
