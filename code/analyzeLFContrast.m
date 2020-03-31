function [] = analyzeLFContrast(subjId)
display(['STARTING - Main Analysis: ',subjId])
% Load the subject relevant info
analysisParams = getSubjectParams(subjId);

analysisParams.preproc = 'hcp';

analysisParams.saveFigs = true;

% Number of bootstrap iterations
numIter  = 200;

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

% Concatenate Packet
[analysisParams, theFullPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);

%% FIT THE TIME COURSE
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theFullPacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCoursePacket = makeDirectionTimeCoursePacketPocket({theFullPacket});

% Make the CRF packet for plotting the stilumi in the nonlinearlty 
stimAndRespForPlot = makeDirectionCrfPacketPocket(analysisParams,iampParams);

% run NR models if set to true
if analysisParams.runNRModels
    % Fit the CRF -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket,'talkToMe',false);
    
    % Fit the CRF with the NR common amplitude -- { } is because this expects a cell
    [~,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket, 'commonAmp', true,'talkToMe',false);
    
    % Fit the CRF with the NR common Exponent -- { } iPs because this expects a cell
    [~,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket, 'commonExp', true,'talkToMe',false);
    
    % Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
    [~,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket, 'commonAmp', true, 'commonExp', true,'talkToMe',false);
end

% Fit the time course with the QCM -- { } is because this expects a cell
[qcmTcOBJ,qcmTcParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000,'talkToMe',false);

%% CRF
% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict CRF from direction model fits
if analysisParams.runNRModels
    % Predict the responses for CRF with params from NR common Amp.
    crfPlot.respNrCrf = nrCrfOBJ.computeResponse(nrCrfParams{1},crfStimulus,[]);
    crfPlot.respNrCrf.plotColor = [.5, .3, .8];
    
    % Predict the responses for CRF with params from NR common Amp.
    crfPlot.respNrCrfAmp = nrCrfOBJ.computeResponse(nrCrfParamsAmp{1},crfStimulus,[]);
    crfPlot.respNrCrfAmp.plotColor = [0, 0, 1];
    
    % Predict the responses for CRF with params from NR common Exp
    crfPlot.respNrCrfExp = nrCrfOBJ.computeResponse(nrCrfParamsExp{1},crfStimulus,[]);
    crfPlot.respNrCrfExp.plotColor = [0, .33, 1];
    
    % Predict the responses for CRF with params from NR common Amp and Exps
    crfPlot.respNrCrfAmpExp = nrCrfOBJ.computeResponse(nrCrfParamsAmpExp{1},crfStimulus,[]);
    crfPlot.respNrCrfAmpExp.plotColor = [0, .66, 1];
end

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmTcOBJ.computeResponse(qcmTcParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.plotColor = qcmColor;

%% TIME COURSE PREDICTIONS

if analysisParams.runNRModels
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    nrAmp = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParamsAmp{1}, timeCoursePacket, 'plotColor', [0, 0, 1]);
    [timeCoursePlot.nrAmp] = chopUpTimeCourse(nrAmp{1},20);
    
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    nrAmpSemi = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParamsExp{1}, timeCoursePacket, 'plotColor', [0, .33, 1]);
    [timeCoursePlot.nrAmpSemi] = chopUpTimeCourse(nrAmpSemi{1},20);
    
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    nrAmpSemiExp = responseFromPacket('nrFullTCPred', analysisParams, nrCrfParamsAmpExp{1}, timeCoursePacket, 'plotColor', [0, .66, 1]);
    [timeCoursePlot.nrAmpSemiExp] = chopUpTimeCourse(nrAmpSemiExp{1},20);
end

% Get the time course predicitions fromt the QCM params fit to the CRF
qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmTcParams{1}, timeCoursePacket, 'plotColor', qcmColor);
[timeCoursePlot.qcm] = chopUpTimeCourse(qcmTimeCourse{1},20);

% Add the IAMP
iampResponses.plotColor = iampColor;
[timeCoursePlot.IAMP]   = chopUpTimeCourse(iampResponses,20);

% Add clean time
theTimeCourse = {theFullPacket.response};
theTimeCourse{1}.plotColor =[0, 0, 0];
[timeCoursePlot.timecourse] = chopUpTimeCourse(theTimeCourse{1},20);

%% Bootstrap
% init vars
iampParamsMat    = [];
iampResponseBoot = [];
qcmParamsMat     = [];
crfQCMBoot       = [];
tcQCMBoot        = [];
for ii = 1:numIter
    [analysisParams, theFullPacketBoot] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',true);
    timeCoursePacketBoot = makeDirectionTimeCoursePacketPocket({theFullPacketBoot});
    
    % fit the IAMP model
    defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
    [iampParamsBoot,fValBoot(ii),iampResps] = iampOBJ.fitResponse(theFullPacketBoot,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    iampParamsMat = [iampParamsMat, iampParamsBoot.paramMainMatrix];
    iampResponseBoot = [iampResponseBoot; iampResps.values];
    
    % Fit the time course with the QCM -- { } is because this expects a cell
    [qcmTcOBJ,qcmCrfParamsBoot] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacketBoot,'fitErrorScalar',1000,'talkToMe',false);
    qcmParamsMat = [qcmParamsMat,qcmTcOBJ.paramsToVec(qcmCrfParamsBoot{1})];
    crfQCMBootStruct = qcmTcOBJ.computeResponse(qcmCrfParamsBoot{1},crfStimulus,[]);
    crfQCMBoot = [crfQCMBoot; crfQCMBootStruct.values];
    
    % get the QCM time course prediction for the bootstrap to calc error bars
    qcmTimeCourseBoot = responseFromPacket('qcmPred', analysisParams, qcmTcParams{1}, timeCoursePacketBoot, 'plotColor', [0, 1, 0]);
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

% Plot the CRF
if analysisParams.showPlots
    
    crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semIAMPParams,...
                     'subtractBaseline', true, 'iampColor',iampColor,'indivBootCRF',[]);
    
    if analysisParams.saveFigs
        figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
            [analysisParams.expSubjID,'_CRF_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        set(crfHndl, 'Renderer', 'Painters');
        FigureSave(figNameCrf,crfHndl,'pdf');
    end
end

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
            [analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        % Save it
        print(tcHndl, figNameTc, '-dpdf', '-r300');
    end
end

%Plot isoresponce contour
if analysisParams.showPlots
 
    eqContrastPts = computeEquivContrast(stimAndRespForPlot,qcmTcParams{1});
    [ellipseNonlinHndl] = plotEllipseAndNonLin(qcmTcParams{1},'plotColor', qcmColor,...
                                 'qcmCI',ciQCMParams,'dispParams',true,'addEqContrastPts', eqContrastPts);
                        
    if analysisParams.saveFigs
        set(ellipseNonlinHndl, 'Renderer', 'Painters');
        figureSizeInches = [13.5 6];
        set(ellipseNonlinHndl, 'PaperUnits', 'inches');
        set(ellipseNonlinHndl, 'PaperSize',figureSizeInches);
        set(ellipseNonlinHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
        figNameEllipseNonlin = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
            [analysisParams.expSubjID,'_Ellipse_Nonlin_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        print(ellipseNonlinHndl, figNameEllipseNonlin, '-dpdf', '-r300');
    end
end
display(['COMPLETED: ',subjId])
end