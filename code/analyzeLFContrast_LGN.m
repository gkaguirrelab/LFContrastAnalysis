function [] = analyzeLFContrast_LGN(subjId)
display(['STARTING - Main Analysis: ',subjId])


%% Loop ever replications
% Load the subject relevant info
numReplications = 2;
for rep = 1:numReplications
    if rep == 2
        subjId = [subjId '_replication'];
    end
    
    analysisParams = getSubjectParams(subjId);
    
    analysisParams.preproc = 'hcp';
    
    % Number of bootstrap iterations
    numIter  = 2;
    
    %turn on or off plotting
    analysisParams.showPlots = true;
    qcmColor  = [0.4078, 0.2784, 0.5765];
    iampColor = [0.8902, 0.6235, 0.5529];
    
    %% Load the relevant data (SDM, HRF, TC)
    
    %set the HRF
    [analysisParams] = loadHRF(analysisParams);
    
    % Load the time course
    [fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);
    
    % Get a packet for each run (1-20)
    [analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData);
    
    % Concatenate Packet
    [analysisParams, theFullPacket{rep}] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);
end
analysisParams.saveFigs = false;
analysisParams.useMedianValidations = false;

%% MAKE THE MEGA PACKET
% make the new timebase
timebase = [theFullPacket{1}.response.timebase,theFullPacket{1}.response.timebase + ...
            (max(theFullPacket{1}.response.timebase)+analysisParams.TR*1000)];
% The response
theMegaPacket.response.values   = [theFullPacket{1}.response.values,theFullPacket{2}.response.values];
theMegaPacket.response.timebase = timebase;
% The stimulus
theMegaPacket.stimulus.values   = [theFullPacket{1}.stimulus.values,theFullPacket{2}.stimulus.values ];
theMegaPacket.stimulus.timebase = timebase;
% The kernel
theMegaPacket.kernel.values     = [((theFullPacket{1}.kernel.values+theFullPacket{2}.kernel.values)./2), ...
                                    zeros(size(theFullPacket{1}.kernel.values))];
theMegaPacket.kernel.timebase   = timebase;
% The metadata
theMegaPacket.metaData.stimDirections = [theFullPacket{1}.metaData.stimDirections,theFullPacket{2}.metaData.stimDirections];
theMegaPacket.metaData.stimContrasts  = [theFullPacket{1}.metaData.stimContrasts,theFullPacket{2}.metaData.stimContrasts];
theMegaPacket.metaData.lmsContrast    = [theFullPacket{1}.metaData.lmsContrast,theFullPacket{2}.metaData.lmsContrast];
 
%% FIT THE TIME COURSE
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(theMegaPacket.stimulus.values,1);
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theMegaPacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCoursePacket = makeDirectionTimeCoursePacketPocket({theMegaPacket});

% Make the CRF packet for plotting the stilumi in the nonlinearlty
stimAndRespForPlot = makeDirectionCrfPacketPocket(analysisParams,iampParams);

% Fit the time course with the QCM -- { } is because this expects a cell
[qcmTcOBJ,qcmTcParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000,'talkToMe',false);

%% CRF
% Upsample the NR repsonses
crfStimulus = upsampleCRF(analysisParams);

% Predict the responses for CRF with params from QCM
crfPlot.respQCMCrf = qcmTcOBJ.computeResponse(qcmTcParams{1},crfStimulus,[]);
crfPlot.respQCMCrf.plotColor = qcmColor;

%% TIME COURSE PREDICTIONS

% Get the time course predicitions fromt the QCM params fit to the CRF
qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmTcParams{1}, timeCoursePacket, 'plotColor', qcmColor);
[timeCoursePlot.qcm] = chopUpTimeCourse(qcmTimeCourse{1},40);

% Add the IAMP
iampResponses.plotColor = iampColor;
[timeCoursePlot.IAMP]   = chopUpTimeCourse(iampResponses,40);

% Add clean time
theTimeCourse = {theMegaPacket.response};
theTimeCourse{1}.plotColor =[0, 0, 0];
[timeCoursePlot.timecourse] = chopUpTimeCourse(theTimeCourse{1},40);

%% Bootstrap
% init vars
iampParamsMat    = [];
iampResponseBoot = [];
qcmParamsMat     = [];
crfQCMBoot       = [];
tcQCMBoot        = [];
for ii = 1:numIter
    [analysisParams, theMegaPacketBoot,runOrder(ii,:)] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',true);
    timeCoursePacketBoot = makeDirectionTimeCoursePacketPocket({theMegaPacketBoot});
    
    % fit the IAMP model
    defaultParamsInfo.nInstances = size(theMegaPacket.stimulus.values,1);
    [iampParamsBoot,fValBoot(ii),~] = iampOBJ.fitResponse(theMegaPacketBoot,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    iampParamsMat = [iampParamsMat, iampParamsBoot.paramMainMatrix];
    
    
    % Fit the time course with the QCM -- { } is because this expects a cell
    [qcmTcOBJ,qcmCrfParamsBoot] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacketBoot,'fitErrorScalar',1000,'talkToMe',false);
    qcmParamsMat = [qcmParamsMat,qcmTcOBJ.paramsToVec(qcmCrfParamsBoot{1})];
    
    % TO PUT ERROR BARS ON THE QCM CRF
    crfQCMBootStruct = qcmTcOBJ.computeResponse(qcmCrfParamsBoot{1},crfStimulus,[]);
    crfQCMBoot = [crfQCMBoot; crfQCMBootStruct.values];
    
    % get the IAMP time course prediction for the bootstrap to calc error bars
    iampResps = computeResponse(iampOBJ,iampParamsBoot,theMegaPacket.stimulus,theMegaPacket.kernel);
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

%contrastType
if analysisParams.showPlots
    
    if analysisParams.useMedianValidations
        figSavePath = fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID,'actualContrast');
    else
        figSavePath =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID);
    end
    
    % Plot the CRF
    crfHndl = plotCRF(analysisParams, crfPlot, crfStimulus, iampParams,semIAMPParams,...
        'subtractBaseline', true, 'iampColor',iampColor,'indivBootCRF',[]);
    
    if analysisParams.saveFigs
        figNameCrf =  fullfile(figSavePath,[analysisParams.expSubjID,'_CRF_' analysisParams.sessionNickname...
            '_' analysisParams.preproc '.pdf']);
        set(crfHndl, 'Renderer', 'Painters');
        FigureSave(figNameCrf,crfHndl,'pdf');
    end
    
    
    % Plot the time course prediction
    tcHndl = plotTimeCourse(analysisParams, timeCoursePlot, zeros(20,1), 20);
    
    if analysisParams.saveFigs
        % Plot configuration
        set(tcHndl, 'Renderer', 'Painters');
        figureSizeInches = [20 15];
        set(tcHndl, 'PaperUnits', 'inches');
        set(tcHndl, 'PaperSize',figureSizeInches);
        set(tcHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
        % Full file name
        figNameTc =  fullfile(figSavePath,[analysisParams.expSubjID,'_TimeCourse_' analysisParams.sessionNickname...
            '_' analysisParams.preproc '.pdf']);
        % Save it
        print(tcHndl, figNameTc, '-dpdf', '-r300');
    end
    
    %Plot isoresponce contour
    eqContrastPts = computeEquivContrast(stimAndRespForPlot,qcmTcParams{1});
    [ellipseNonlinHndl] = plotEllipseAndNonLin(qcmTcParams{1},'plotColor', qcmColor,...
        'qcmCI',ciQCMParams,'dispParams',true,'addEqContrastPts', eqContrastPts);
    
    if analysisParams.saveFigs
        set(ellipseNonlinHndl, 'Renderer', 'Painters');
        figureSizeInches = [13.5 6];
        set(ellipseNonlinHndl, 'PaperUnits', 'inches');
        set(ellipseNonlinHndl, 'PaperSize',figureSizeInches);
        set(ellipseNonlinHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
        figNameEllipseNonlin = fullfile(figSavePath,[analysisParams.expSubjID,'_Ellipse_Nonlin_' ...
            analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
        print(ellipseNonlinHndl, figNameEllipseNonlin, '-dpdf', '-r300');
    end
end

%% Check the CI on the timecourse -- Ploting run 1
figure; hold on
plot(timeCoursePlot.timecourse{1}.timebase, [tcQCMBoot(:,1:360)]','Color',[.1 .4 .8 .5])
plot(timeCoursePlot.qcm{1}.timebase,timeCoursePlot.qcm{1}.values, 'Color',[.1 .2 1],'LineWidth',2)
display(['COMPLETED: ',subjId])
end