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


%% TIME COURSE PREDICTIONS

% Get the time course predicitions fromt the QCM params fit to the CRF
qcmPredSessionOne = responseFromPacket('qcmPred', analysisParams, qcmParamsSessionTwo{1}, timeCoursePacketSessionOne, 'plotColor', qcmColor);
qcmPredSessionTwo = responseFromPacket('qcmPred', analysisParams, qcmParamsSessionOne{1}, timeCoursePacketSessionTwo, 'plotColor', qcmColor);

qcmtimeCourseSessionOne = chopUpTimeCourse(qcmPredSessionOne{1},10);
qcmtimeCourseSessionTwo = chopUpTimeCourse(qcmPredSessionTwo{1},10);
timeCoursePlot.qcm   = [qcmtimeCourseSessionOne,qcmtimeCourseSessionTwo];
% Add clean time
theTimeCourseSessionOne = {sessionOnePacket.response};
theTimeCourseSessionOne{1}.plotColor =[0, 0, 0];

theTimeCourseSessionTwo = {sessionTwoPacket.response};
theTimeCourseSessionTwo{1}.plotColor =[0, 0, 0];

timeCoursePlotSessionOne = chopUpTimeCourse(theTimeCourseSessionOne{1},10);
timeCoursePlotSessionTwo = chopUpTimeCourse(theTimeCourseSessionTwo{1},10);
timeCoursePlot.timecourse = [timeCoursePlotSessionOne,timeCoursePlotSessionTwo];

%% MAKE THE PLOTS

% Plot the time course prediction
if analysisParams.showPlots
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
            [analysisParams.expSubjID,'_SplitSession_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        % Save it
        print(tcHndl, figNameTc, '-dpdf', '-r300');
    end
end

