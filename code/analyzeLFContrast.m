%% Set up params
% Set the subject: 'LZ23', 'KAS25', 'AP26'
subjId = 'KAS25';

% Number of bootstrap iterations
numIter  = 100;

% Load the subject relevant info
analysisParams = getSubjectParams(subjId);

% Flag for running all the NR models
analysisParams.runNRModels = false;

%set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

%turn on or off plotting
analysisParams.showPlots = true;

%% Load the relevant data (SDM, HRF, TC)

%set the HRF
[analysisParams] = loadHRF(analysisParams);

% Load the time course
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

%% Get a packet for each run (1-20)
[analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData);

%% Fit the GLM to get the crf betas
[analysisParams, theFullPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);

% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theFullPacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

%% Bootstrap
% init vars
paramsMat   = [];
for ii = 1:numIter
    [analysisParams, theFullPacketBoot] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',true);
    
    % fit the IAMP model
    defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
    [iampParamsBoot,fValBoot(ii),~] = iampOBJ.fitResponse(theFullPacketBoot,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    paramsMat = [paramsMat, iampParamsBoot.paramMainMatrix];
end

% Calc SEM for betas
semIAMP = std(paramsMat,0,2);
semParams = iampParams;
semParams.paramMainMatrix =semIAMP;

%% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCoursePacket = makeDirectionTimeCoursePacketPocket({theFullPacket});

% Make the CRF packet
directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampParams);

%% Fit the direction based models to the concatenated time course

% run NR models if set to true
if analysisParams.runNRModels
    % Fit the CRF -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket);
    
    % Fit the CRF with the NR common amplitude -- { } is because this expects a cell
    [~,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket, 'commonAmp', true);
    
    % Fit the CRF with the NR common Exponent -- { } iPs because this expects a cell
    [~,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket, 'commonExp', true);
    
    % Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
    [~,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', timeCoursePacket, 'commonAmp', true, 'commonExp', true);
end

% Fit the time course with the QCM -- { } is because this expects a cell
[qcmCrfOBJ,qcmCrfParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000);

% Do some plotting of these fits
if analysisParams.showPlots
    [nrVals] = plotNakaRushtonFromParams(qcmCrfParams{1}.crfAmp ,qcmCrfParams{1}.crfExponent,qcmCrfParams{1}.crfSemi,...
        'analysisParams',analysisParams,'plotFunction',true,'savePlot',true);
end

% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict CRF from direction model fits
if analysisParams.runNRModels
    % Predict the responses for CRF with params from NR common Amp.
    crfPlot.respNrCrf = nrCrfOBJ.computeResponse(nrCrfParams{1},crfStimulus,[]);
    crfPlot.respNrCrf.color = [.5, .3, .8];
    
    % Predict the responses for CRF with params from NR common Amp.
    crfPlot.respNrCrfAmp = nrCrfOBJ.computeResponse(nrCrfParamsAmp{1},crfStimulus,[]);
    crfPlot.respNrCrfAmp.color = [0, 0, 1];
    
    % Predict the responses for CRF with params from NR common Exp
    crfPlot.respNrCrfExp = nrCrfOBJ.computeResponse(nrCrfParamsExp{1},crfStimulus,[]);
    crfPlot.respNrCrfExp.color = [0, .33, 1];
    
    % Predict the responses for CRF with params from NR common Amp and Exps
    crfPlot.respNrCrfAmpExp = nrCrfOBJ.computeResponse(nrCrfParamsAmpExp{1},crfStimulus,[]);
    crfPlot.respNrCrfAmpExp.color = [0, .66, 1];
end

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmCrfOBJ.computeResponse(qcmCrfParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.color = [0, 1, 0];

% Plot the CRF from the IAMP, QCM, and  fits
if analysisParams.showPlots
    crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semParams,'subtractBaseline', true);
    figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_CRF_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    FigureSave(figNameCrf,crfHndl,'pdf');
end

% Get the time course predicitions of the CRF params

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
qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmCrfParams{1}, timeCoursePacket, 'plotColor', [0, 1, 0]);
[timeCoursePlot.qcm] = chopUpTimeCourse(qcmTimeCourse{1},20);

% Add the IAMP
iampResponses.plotColor = [.5,.5,.5];
[timeCoursePlot.IAMP] = chopUpTimeCourse(iampResponses,20);

% Add clean time
theTimeCourse = {theFullPacket.response};
theTimeCourse{1}.plotColor =[0, 0, 0];
[timeCoursePlot.timecourse] = chopUpTimeCourse(theTimeCourse{1},20);

% Plot the time course prediction for each run using the different fits to
% the crf
if analysisParams.showPlots
    tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, zeros(20,1), 20);
    figNameTc =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    FigureSave(figNameTc,tcHndl,'pdf');
end

% Plot isoresponce contour
% if analysisParams.showPlots
%     thresholds = [0.15];
%     colors     = [0.5,0.0,0.0];
%     qcmHndl    = plotIsoresponse(analysisParams,iampParams,qcmCrfMeanParams,thresholds,nrCrfParamsAmp,colors);
% %     figNameQcm = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
% %         [analysisParams.expSubjID,'_QCM_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
% %     FigureSave(figNameQcm,qcmHndl,'pdf');
% end
