function [] = analyzeLFContrast(subjId)
display(['STARTING - Main Analysis: ',subjId])
% Load the subject relevant info
analysisParams = getSubjectParams(subjId);

analysisParams.preproc = 'hcp';

analysisParams.saveFigs = true;

% Flag for running all the NR models
analysisParams.runNRModels = false;

% bandpass the signal
analysisParams.highpass = false;


% number of bootstrap iterations
numIter = 200;

%turn on or off plotting
analysisParams.showPlots = true;
qcmColor  = [0.4078, 0.2784, 0.5765];
iampColor = [0.8902, 0.6235, 0.5529];
splitSessionColor = [.27 .58 .12];

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

%% FULL SESSION IAMP
% make the full session packet
[analysisParams, theFullPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);
%% FIT THE TIME COURSE
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theFullPacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

% Fit the time course with the QCM -- { } is because this expects a cell
timeCoursePacket = makeDirectionTimeCoursePacketPocket({theFullPacket});
[qcmTcOBJ,qcmTcParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000,'talkToMe',false);
crfStimulus = upsampleCRF(analysisParams);
%dummy up sem for iamp to be 0
semIAMPParams = iampParams;
semIAMPParams.paramMainMatrix = zeros(size(iampParams.paramMainMatrix))';

crfPlot.QcmCrf = qcmTcOBJ.computeResponse(qcmTcParams{1},crfStimulus,[]);
crfPlot.QcmCrf.plotColor = qcmColor;
%% SPLIT SESSION

% Split the sessions and concatenate Packet
[analysisParams, sessionOnePacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket(1:10),'bootstrap',false);
[analysisParams, sessionTwoPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket(11:20),'bootstrap',false);
%% FIT THE TIME COURSE WITH THE QCM

% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCoursePacketSessionOne = makeDirectionTimeCoursePacketPocket({sessionOnePacket});
timeCoursePacketSessionTwo = makeDirectionTimeCoursePacketPocket({sessionTwoPacket});

% Fit the time course with the QCM -- { } is because this expects a cell
[qcmTcOBJ,qcmParamsSessionOne] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacketSessionOne,'fitErrorScalar',1000,'talkToMe',false);
[qcmTcOBJ,qcmParamsSessionTwo] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacketSessionTwo,'fitErrorScalar',1000,'talkToMe',false);

%% Contrast Resposne functions
% Upsample the NR repsonses

crfStimulusSessionOne.values = crfStimulus.values(:,1:100);
crfStimulusSessionOne.timebase = crfStimulus.timebase(1:100);
crfStimulusSessionTwo.values = crfStimulus.values(:,101:200);
crfStimulusSessionTwo.timebase = crfStimulus.timebase(101:200);


% Predict the responses for CRF with params from QCM
crfRespQcmSessionOne = qcmTcOBJ.computeResponse(qcmParamsSessionTwo{1},crfStimulusSessionOne,[]);
crfRespQcmSessionTwo = qcmTcOBJ.computeResponse(qcmParamsSessionOne{1},crfStimulusSessionTwo,[]);

crfPlot.QcmCrfSplit.values = [crfRespQcmSessionOne.values,crfRespQcmSessionTwo.values];
crfPlot.QcmCrfSplit.timebase = [crfRespQcmSessionOne.timebase,crfRespQcmSessionTwo.timebase];
crfPlot.QcmCrfSplit.plotColor = splitSessionColor;

%% TIME COURSE PREDICTIONS
% Get the time course predicitions fromt the QCM params fit to the CRF
qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmTcParams{1}, timeCoursePacket, 'plotColor', qcmColor);
[timeCoursePlot.qcm] = chopUpTimeCourse(qcmTimeCourse{1},20);

% Get the time course predicitions fromt the QCM params fit to the CRF
qcmPredSessionOne = responseFromPacket('qcmPred', analysisParams, qcmParamsSessionTwo{1}, timeCoursePacketSessionOne, 'plotColor', splitSessionColor);
qcmPredSessionTwo = responseFromPacket('qcmPred', analysisParams, qcmParamsSessionOne{1}, timeCoursePacketSessionTwo, 'plotColor', splitSessionColor);

qcmtimeCourseSessionOne = chopUpTimeCourse(qcmPredSessionOne{1},10);
qcmtimeCourseSessionTwo = chopUpTimeCourse(qcmPredSessionTwo{1},10);
timeCoursePlot.qcmSplit   = [qcmtimeCourseSessionOne,qcmtimeCourseSessionTwo];
% Add clean time
theTimeCourseSessionOne = {sessionOnePacket.response};
theTimeCourseSessionOne{1}.plotColor =[0, 0, 0];

theTimeCourseSessionTwo = {sessionTwoPacket.response};
theTimeCourseSessionTwo{1}.plotColor =[0, 0, 0];

timeCoursePlotSessionOne = chopUpTimeCourse(theTimeCourseSessionOne{1},10);
timeCoursePlotSessionTwo = chopUpTimeCourse(theTimeCourseSessionTwo{1},10);
timeCoursePlot.timecourse = [timeCoursePlotSessionOne,timeCoursePlotSessionTwo];


%% Bootstrap
% init vars
iampParamsMat    = [];
iampResponseBoot = [];
qcmParamsMat     = [];
crfQCMBoot       = [];
tcQCMBoot        = [];
for ii = 1:numIter
    % Session One
    [analysisParams, theFullPacketBootOne,runOrderOne(ii,:)] = concatPackets(analysisParams, iampTimeCoursePacketPocket{1:10},'bootstrap',true);
    timeCoursePacketBootOne = makeDirectionTimeCoursePacketPocket({theFullPacketBoot});
    
    % Session Two
    [analysisParams, theFullPacketBootTwo,runOrderTwo(ii,:)] = concatPackets(analysisParams, iampTimeCoursePacketPocket{11:20},'bootstrap',true);
    timeCoursePacketBootTwo = makeDirectionTimeCoursePacketPocket({theFullPacketBoot});
    
    
    % Fit the time course with the QCM session one packets
    [qcmTcOBJ,qcmCrfParamsBoot] = fitDirectionModel(analysisParams, 'qcmFit', theFullPacketBootOne,'fitErrorScalar',1000,'talkToMe',false);
    qcmParamsMat = [qcmParamsMat,qcmTcOBJ.paramsToVec(qcmCrfParamsBoot{1})];
    
    
    % Fit the time course with the QCM session two packets
    [qcmTcOBJ,qcmCrfParamsBoot] = fitDirectionModel(analysisParams, 'qcmFit', theFullPacketBootTwo,'fitErrorScalar',1000,'talkToMe',false);
    qcmParamsMat = [qcmParamsMat,qcmTcOBJ.paramsToVec(qcmCrfParamsBoot{1})];
    
    % TO PUT ERROR BARS ON THE QCM CRF
    crfQCMBootStruct = qcmTcOBJ.computeResponse(qcmCrfParamsBoot{1},crfStimulus,[]);
    crfQCMBoot = [crfQCMBoot; crfQCMBootStruct.values];
    
    % get the IAMP time course prediction for the bootstrap to calc error bars
    iampResps = computeResponse(iampOBJ,iampParamsBoot,theFullPacket.stimulus,theFullPacket.kernel);
    iampResponseBoot = [iampResponseBoot; iampResps.values];
    
    % get the QCM time course prediction for the bootstrap to calc error bars
    qcmTimeCourseBoot = responseFromPacket('qcmPred', analysisParams, qcmCrfParamsBoot{1}, timeCoursePacket, 'plotColor', [0, 1, 0]);
    tcQCMBoot         = [tcQCMBoot; qcmTimeCourseBoot{1}.values];
end

%% Get error bars on stuff
ciPercent = 68;
upperCiVal = 100 - ((100 - ciPercent)/2);
lowerCiVal = ((100 - ciPercent)/2);

% Calc Error for IAMP betas
ciIampParams = prctile(iampParamsMat',[upperCiVal lowerCiVal]);
ciIampParams = abs(ciIampParams - iampParams.paramMainMatrix');
semIAMPParams = iampParams;
semIAMPParams.paramMainMatrix =ciIampParams;

% Calc error for IAMP TC Prediction
ciIampTC = prctile(iampResponseBoot,[upperCiVal lowerCiVal]);
errorIampTC = abs(ciIampTC -qcmTimeCourse{1}.values);
timeCoursePlot.IAMP = addErrorBarsToTimeCouse(errorIampTC,timeCoursePlot.IAMP);

% Calc error for QCM Params
ciQCM = prctile(qcmParamsMat',[upperCiVal lowerCiVal])';
ciQCMParams.mar   =ciQCM(1,:);
ciQCMParams.angle =ciQCM(2,:);
ciQCMParams.amp   =ciQCM(3,:);
ciQCMParams.exp   =ciQCM(4,:);
ciQCMParams.semi  =ciQCM(5,:);

% Calc error for QCM CRF
respQCMCrfCI = prctile(crfQCMBoot,[upperCiVal lowerCiVal]);
crfPlot.respQCMCrf.shaddedErrorBars  = abs(respQCMCrfCI -crfPlot.respQCMCrf.values);

% Calc SEM for QCM TC Prediction
ciQcmTC =prctile(tcQCMBoot,[upperCiVal lowerCiVal]);
errorQcmTC = abs(ciQcmTC -qcmTimeCourse{1}.values);
timeCoursePlot.qcm = addErrorBarsToTimeCouse(errorQcmTC,timeCoursePlot.qcm);


%% MAKE THE PLOTS

% Plot the time course prediction
if analysisParams.showPlots
    crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semIAMPParams,...
        'subtractBaseline', true, 'iampColor',iampColor,'indivBootCRF',[]);
    
    if analysisParams.saveFigs
        figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
            [analysisParams.expSubjID,'_SplitSession_CRF_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        set(crfHndl, 'Renderer', 'Painters');
        FigureSave(figNameCrf,crfHndl,'pdf');
    end
    
    tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, zeros(20,1), 20);
    
    if analysisParams.saveFigs
        % Plot configuration
        set(tcHndl, 'Renderer', 'Painters');
        figureSizeInches = [20 15];
        set(tcHndl, 'PaperUnits', 'inches');
        set(tcHndl, 'PaperSize',figureSizeInches);
        set(tcHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
        % Full file name
        figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
            [analysisParams.expSubjID,'_SplitSession_TC_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        % Save it
        print(tcHndl, figNameTc, '-dpdf', '-r300');
    end
    
end

