%crossValidateDirectionModels
heldOutRunOrder = [5, 10, 3, 7, 2, 8, 6, 9, 4, 1; ...
    10, 3, 4, 8, 1, 6, 7, 9, 2, 5];

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('LZ23');

% Clip fisrt 2 TRs from time series?
% if no clipping then put 0;
analysisParams.numClipFramesStart = 0;
analysisParams.numClipFramesEnd   = 2;

% Make mask from the area and eccentricity maps
analysisParams.areaNum     = 1;
analysisParams.eccenRange  = [1 20];

% Define the TR
analysisParams.TR = 0.800;
analysisParams.baselineCondNum = 6;
analysisParams.timeStep = 1/100;
analysisParams.generateIAMPPlots = false;
analysisParams.generateCrossValPlots = false;

% Plotting params
analysisParams.numSamples = 25;

% Get the cleaned time series
[fullCleanData, analysisParams] = getTimeCourse(analysisParams);

[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% create the random draws with replacement
sampleMatrix = randi([1,10],length(analysisParams.sessionFolderName),analysisParams.numAcquisitions);



% Cross validation folds
for ii = 1:size(heldOutRunOrder,2)
    
    heldOutParams = [iampParams(1,heldOutRunOrder(1,ii)); iampParams(2,heldOutRunOrder(2,ii))];
    heldOutTCPackets = [iampTimeCoursePacketPocket(1,heldOutRunOrder(1,ii)) ;iampTimeCoursePacketPocket(2,heldOutRunOrder(2,ii))];
    heldOutRawTC = [rawTC(1,heldOutRunOrder(1,ii)); rawTC(2,heldOutRunOrder(2,ii))];
    heldOutBaselineShift = [concatBaselineShift(1,heldOutRunOrder(1,ii)); concatBaselineShift(2,heldOutRunOrder(2,ii))];
    
    % get held out params -- i am sure this is a bad way to code this but
    % it is late
    tmpMat = ones(size(iampParams));
    tmpMat(1,heldOutRunOrder(1,ii)) = 0;
    tmpMat(2,heldOutRunOrder(2,ii)) = 0;
    cvParams(1,:) = iampParams(logical(tmpMat(1,:)));
    cvParams(2,:) = iampParams(logical(tmpMat(2,:)));
    
    cvTCPackets(1,:) = iampTimeCoursePacketPocket(logical(tmpMat(1,:)));
    cvTCPackets(2,:) = iampTimeCoursePacketPocket(logical(tmpMat(2,:)));
    
    % get fits
    % Get directon/contrast form of time course and IAMP crf packet pockets.
    %
    % This conversion is possible because the IAMP packet pocket has meta data
    % that we put there to allow exactly this conversion.  That meta data
    % encapsulates the key things we need to know about the stimulus obtained
    % from the analysis parameters.
    directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(heldOutTCPackets);
    
    % This puts together pairs of acquistions from the two sessions, so that
    % we have one IAMP fit for each pair.  We do this because to fit the
    % quadratic model, we need data for all of the color directions together.
    %
    % NOTE: This bit is very specific to the design of the experiment we are
    % currently analyzing, and has to do specifically with the way color
    % directions were studied across acquisitions and sessions.
    for jj = 1:analysisParams.numAcquisitions-1
        [concatParams{ii},concatBaselineShift(:,jj)] = iampOBJ.concatenateParams(cvParams(:,jj),'baselineMethod','makeBaselineZero');
    end
    
    directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampOBJ.averageParams(concatParams));
    
    %% Fit the direction based models to the mean IAMP beta weights
    %
    % Fit the CRF with the NR -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket});
    
    % Fit the CRF with the NR common amplitude -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsAmp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true);
    
    % Fit the CRF with the NR common Exponent -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonExp', true);
    
    % Fit the CRF with the NR common amplitude, and exponent  -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsAmpExp] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonExp', true);
    
    % Fit the CRF with the QCM -- { } is because this expects a cell
    [qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});
    
    
    %% Get the time course predicitions of the CRF params
    %
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    timeCoursePlot.nr = responseFromPacket('nrPred', analysisParams, nrCrfParams{1}, directionTimeCoursePacketPocket, 'plotColor', [0, 0, 1]);
    
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    timeCoursePlot.nrAmp = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmp{1}, directionTimeCoursePacketPocket, 'plotColor', [.3, .2, .8]);
    
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    timeCoursePlot.nrExp = responseFromPacket('nrPred', analysisParams, nrCrfParamsExp{1}, directionTimeCoursePacketPocket, 'plotColor', [.5, .66, .9]);
    
    % Get the time course predicitions from the NR common Amp and Semi fit to the CRF
    timeCoursePlot.nrAmpExp = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmpExp{1}, directionTimeCoursePacketPocket, 'plotColor', [0.7, .6, .1]);
    
    % Get the time course predicitions fromt the QCM params fit to the CRF
    timeCoursePlot.qcm = responseFromPacket('qcmPred', analysisParams, qcmCrfMeanParams{1}, directionTimeCoursePacketPocket, 'plotColor', [0, 1, 0]);
    
    % Add clean time
    timeCoursePlot.heldOutRawTC = heldOutRawTC;
    
    % %Plot the time course prediction for each run using the different fits to
    % %the crf
    %plotTimeCourse(analysisParams, timeCoursePlot, heldOutBaselineShift);
    
    %% Calculate RMSE
    for kk = 1:length(timeCoursePlot.heldOutRawTC)
        
        % Calculate the RMSE for the  NR preds
        nrRMSE(kk,ii) = sqrt(mean((timeCoursePlot.nr{kk}.values-timeCoursePlot.heldOutRawTC{1}.values).^2));
        
        % Calculate the RMSE for the NR common Amp preds
        nrAmpRMSE(kk,ii) = sqrt(mean((timeCoursePlot.nrAmp{kk}.values-timeCoursePlot.heldOutRawTC{1}.values).^2));
        
        % Calculate the RMSE for the NR common Exp preds
        nrExpRMSE(kk,ii) = sqrt(mean((timeCoursePlot.nrExp{kk}.values-timeCoursePlot.heldOutRawTC{1}.values).^2));
        
        % Calculate the RMSE for the NR common Amp and Exp preds
        nrAmpExpRMSE(kk,ii) = sqrt(mean((timeCoursePlot.nrAmpExp{kk}.values-timeCoursePlot.heldOutRawTC{1}.values).^2));
        
        % Calculate the RMSE for the qcm
        qcmRMSE(kk,ii) = sqrt(mean((timeCoursePlot.qcm{kk}.values-timeCoursePlot.heldOutRawTC{1}.values).^2));
        
    end
    
    
end



%% Get the mean and SEM of the crossval

meanNrRMSE = mean(nrRMSE(:));
semNrRMSE = std(nrRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common Amp preds
meanNrAmpRMSE = mean(nrAmpRMSE(:));
semNrAmpRMSE= std(nrAmpRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common Exp preds
meanNrExpRMSE = mean(nrExpRMSE(:));
semNrExpRMSE = std(nrExpRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common Amp and Exp preds
meanNrAmpExpRMSE = mean(nrAmpExpRMSE(:));
semNrAmpExpRMSE = std(nrAmpExpRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the qcm
meanQcmRMSE= mean(qcmRMSE(:));
semQcmRMSE = std(qcmRMSE(:))./sqrt(size(heldOutRunOrder,2));

plotnames = categorical({'meanNrRMSE','meanNrAmpRMSE','meanNrExpRMSE','meanNrAmpExpRMSE','meanQcmRMSE'});
barsForPlot = [meanNrRMSE,meanNrAmpRMSE,meanNrExpRMSE,meanNrAmpExpRMSE,meanQcmRMSE];
errorForPlot = [semNrRMSE,semNrAmpRMSE,semNrExpRMSE,semNrAmpExpRMSE,semQcmRMSE];
bar(barsForPlot)
errorbar(barsForPlot,errorForPlot)
