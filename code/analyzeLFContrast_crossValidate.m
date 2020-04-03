function [] = analyzeLFContrast_crossValidate(subjId)

display(['STARTING - Cross Validating: ',subjId])
% Load the subject relevant info
analysisParams = getSubjectParams(subjId);

analysisParams.preproc = 'hcp';

analysisParams.saveFigs = true;

% bandpass the signal
analysisParams.highpass = false;

%turn on or off plotting
analysisParams.showPlots = true;
qcmColor  = [0.4078, 0.2784, 0.5765];
iampColor = [0.8902, 0.6235, 0.5529];

%% Load the relevant data (SDM, HRF, TC)

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

%% unscramble the order
ordr = leaveOutPairs';
ordr = ordr(:)';
[~,indx] = sort(ordr);


timeCoursePlot.qcm = timeCoursePlot.qcm(indx);
timeCoursePlot.IAMP = timeCoursePlot.IAMP(indx);
timeCoursePlot.timecourse = timeCoursePlot.timecourse(indx);
rSquaredQcm = rSquaredQcm(indx);
rSquaredIamp = rSquaredIamp(indx);

%% Calc mean and CI for R squared
rSquaredQcmMean = mean(rSquaredQcm);
rSquaredIampMean = mean(rSquaredIamp);

ciPercent = 68;
upperCiVal = 100 - ((100 - ciPercent)/2);
lowerCiVal = ((100 - ciPercent)/2);

qcmRsquaredCI =prctile(rSquaredQcm,[upperCiVal lowerCiVal]);
iampRsquaredCI =prctile(rSquaredIamp,[upperCiVal lowerCiVal]);

%% MAKE THE PLOTS
tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, zeros(20,1), 20);
if analysisParams.saveFigs
    set(tcHndl, 'Renderer', 'Painters');
    figureSizeInches = [20 15];
    set(tcHndl, 'PaperUnits', 'inches');
    set(tcHndl, 'PaperSize',figureSizeInches);
    set(tcHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
    figNameEllipseNonlin = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_TimeCourse_Cross_Val_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    print(tcHndl, figNameEllipseNonlin, '-dpdf', '-r300');
end

%% Plot the
crossValR2 = figure; hold on;
set(gca,'Box', 'off','linewidth',3,'FontSize',12);
X = categorical({'QCM','GLM'});
X = reordercats(X,{'QCM','GLM'});
b = bar(X,[rSquaredQcmMean;rSquaredIampMean]);
er = errorbar(X,[rSquaredQcmMean;rSquaredIampMean],[qcmRsquaredCI(2),iampRsquaredCI(2)], [qcmRsquaredCI(1),iampRsquaredCI(1)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;
b.FaceColor = 'flat';
b.CData(1,:) = iampColor;
b.CData(2,:) = qcmColor;
b.EdgeColor = [0,0,0];
b.LineWidth = 2;
hXLabel = xlabel('Models');
hYLabel = ylabel('R^{2}');
hTitle  = title('Cross Validated R^{2}');
set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 14);
set( hTitle, 'FontSize', 14,'FontWeight' , 'bold');
ylim([0,1])
set(gca,'TickDir', 'out');
set(gcf, 'Color', 'white' );

modelTxtQcm   = ['{QCM: ' num2str(round(rSquaredQcmMean,2))...
    ' CI [' num2str(round(qcmRsquaredCI(2),2)) ', ' num2str(round(qcmRsquaredCI(1),2)) ']}'];
modelTxtIamp  = ['{IAMP: ' num2str(round(rSquaredIampMean,2))...
    ' CI [' num2str(round(iampRsquaredCI(2),2)) ', ' num2str(round(iampRsquaredCI(1),2)) ']}'];
theTextHandle = text(gca, .55,.95 , modelTxtQcm, 'Interpreter', 'latex');
set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);
theTextHandle = text(gca, .55,0.87, modelTxtIamp, 'Interpreter', 'latex');
set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);


if analysisParams.saveFigs
    set(crossValR2, 'Renderer', 'Painters');
    figureSizeInches = [6 5.5];
    set(crossValR2, 'PaperUnits', 'inches');
    set(crossValR2, 'PaperSize',figureSizeInches);
    set(crossValR2, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
    figNameEllipseNonlin = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_Cross_Val_R2_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    print(crossValR2, figNameEllipseNonlin, '-dpdf', '-r300');
end
display(['COMPLETED: ',subjId])
end
