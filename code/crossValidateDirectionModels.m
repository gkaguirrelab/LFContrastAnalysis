% Cross validate models for LFContrast using both the time course and the crf. 
%
% Description:
%   This script runs the main cross-validated analysis for the QCM and
%   related models on our 2018-19 fMRI color contrast response function
%   datasets.

% History:
%   05/12/19  dhb  More comments and some variable renaming.

%% Clear and close
clear; close all;

%% Define which subject we'll run
%
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('LZ23');
   
%% Set up cross-validation order
%
% For each subject, we have 10 repeats for each stimulus
% set, and need to combine in pairs to have a whole dataset
% held out.  We simply fix that, and specify it here.
heldOutRunOrder = [5, 10, 3, 7, 2, 8, 6, 9, 4, 1; ...
                   10, 3, 4, 8, 1, 6, 7, 9, 2, 5];

%% Set key analysis parameters and store in structure
%
% Clip first/last TRs from time series?
% If no clipping then put 0;
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

%% Get the cleaned time series.
%
% These have been pre-processed and stored for analysis here.
[fullCleanData, analysisParams] = getTimeCourse(analysisParams);

%% Fit the general linear model (IAMP model) to the time series
%
% We fit the IAMP to each run. In the returned arrays, the rows correspond
% to the two different stimulus sets, and the columns to the individual
% runs.
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

%% Handle fact that baseline varies across sessions.
%
% The mean BOLD response varies across sessions and is modeled by the IAMP.
% Here we find the amount each session's data gets shifted by to bring the
% baseline response (response to no stimulus) to zero.  This allows us
% to concatenate data across sessions in a sensible manner below.
for ii = 1:analysisParams.numAcquisitions
    [~,concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end

%% Cross validation folds
%
% In each fold, put together data from N-1 runs, fit, and
% evaluate fit against the held out run.
for ii = 1:size(heldOutRunOrder,2)
    
    %% Get the parameters and data from the held out run
    heldOutParams = [iampParams(1,heldOutRunOrder(1,ii)) ; iampParams(2,heldOutRunOrder(2,ii))];
    heldOutTCPackets = [iampTimeCoursePacketPocket(1,heldOutRunOrder(1,ii)) ; iampTimeCoursePacketPocket(2,heldOutRunOrder(2,ii))];
    heldOutRawTC = [rawTC(1,heldOutRunOrder(1,ii)) ; rawTC(2,heldOutRunOrder(2,ii))];
    heldOutBaselineShift = [concatBaselineShift(1,heldOutRunOrder(1,ii)) ; concatBaselineShift(2,heldOutRunOrder(2,ii))];
    
    %% Get parameters of runs that will be fit in this fold (we'll call
    % these the CV runs).
    %
    % See setdiff() for a possible slicker way to do this.
    %
    % Get a matrix of ones and set entries corresponding to held out runs
    % to zero.
    indxMat = ones(size(iampParams));
    indxMat(1,heldOutRunOrder(1,ii)) = 0;
    indxMat(2,heldOutRunOrder(2,ii)) = 0;
    
    %% Pull out the IAMP parameters and packet pockets for the runs being
    % fit in this fold.
    cvIampParams(1,:) = iampParams(1,logical(indxMat(1,:)));
    cvIampParams(2,:) = iampParams(2,logical(indxMat(2,:)));
    cvIampTimeCoursePacketPockets(1,:) = iampTimeCoursePacketPocket(1,logical(indxMat(1,:)));
    cvIampTimeCoursePacketPockets(2,:) = iampTimeCoursePacketPocket(2,logical(indxMat(2,:)));
    
    %% Create the mean time course model from the cv runs. This just fits
    % all the runs with a constant (over time) response whose value is the
    % mean value across all the runs.  We had better do better than this
    % with our real models.
    for pp = 1:size(cvIampTimeCoursePacketPockets,1)
        for mm = 1:size(cvIampTimeCoursePacketPockets,2)
            cvMeanResponseVal(pp,mm) =  mean(cvIampTimeCoursePacketPockets{pp,mm}.response.values);
        end
    end
    cvMeanModel = ones(size(cvIampTimeCoursePacketPockets{1}.response.values)).*mean(cvMeanResponseVal(:));
    
    %% Get fits for models of interest to the cv data.
    %
    % First together pairs of acquistions from the two sessions, so that
    % we have one IAMP fit for each pair.  We do this because to fit the
    % quadratic model, we need data for all of the color directions together.
    %
    % NOTE: This bit is very specific to the design of the experiment we are
    % currently analyzing, and has to do specifically with the way color
    % directions were studied across acquisitions and sessions.
    for jj = 1:analysisParams.numAcquisitions-1
        [concatCvIampParams{jj},~] = iampOBJ.concatenateParams(cvIampParams(:,jj),'baselineMethod','makeBaselineZero');
    end
    
    % Make the mean IAMP CRF packet by averaging the IAMP parameters over
    % the CV runs.
    directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,iampOBJ.averageParams(concatCvIampParams));
    
    %% Fit the direction based models to the mean IAMP beta weights
    %
    % Fit the CRF with the NR -- { } is because this expects a cell. This
    % has common offset.
    [nrCrfOBJ,nrCrfParams,nrObjFitResponses] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket});
    
    % Fit the CRF with the NR common amplitude -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsAmp, nrAmpObjFitResponses] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true);
    
    % Fit the CRF with the NR common exponent -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsExp,nrExpObjFitResponses] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonExp', true);
    
    % Fit the CRF with the NR common amplitude, and exponent -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsAmpExp, nrAmpExpObjFitResponses] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonExp', true);
    
    % Fit the CRF with the NR common amplitude, exponent,  and semi -- { } is because this expects a cell
    [nrCrfOBJ,nrCrfParamsAmpExpSemi,nrAmpExpSemiObjFitResponses] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket}, 'commonAmp', true, 'commonExp', true, 'commonSemi', true);
    
    % Fit the CRF with the QCM -- { } is because this expects a cell
    [qcmCrfMeanOBJ,qcmCrfMeanParams, qcmObjFitResponses] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});
    
    %% Calculate held out RMSE for the CRF, given CV fits above
    %
    % Concat IAMP params for held out data so we can compare to this.
    [concatHeldOutIampParams{ii},~] = iampOBJ.concatenateParams(heldOutParams,'baselineMethod','makeBaselineZero');
    
    % Mean model 
    meanCrfResp = ones(size(directionCrfMeanPacket.response.values)).*mean(directionCrfMeanPacket.response.values);
    mmCrfRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-meanCrfResp).^2));
    
    % Naka-Rushton with common offset
    nrCrfRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-nrObjFitResponses{1}.values).^2));
    
    % Naka-Rushton with common amp
    nrCrfAmpRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-nrAmpObjFitResponses{1}.values).^2));
    
    % Naka-Rushton with common exp
    nrCrfExpRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-nrExpObjFitResponses{1}.values).^2));
    
    % Naka-Rushton with common amp, exp
    nrCrfAmpExpRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-nrAmpExpObjFitResponses{1}.values).^2));
   
    % Naka-Rushton with common amp, exp, semi
    nrCrfAmpExpSemiRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-nrAmpExpSemiObjFitResponses{1}.values).^2));
    
    % QCM 
    qcmCrfRMSE(ii) = sqrt(mean((concatHeldOutIampParams{ii}.paramMainMatrix-qcmObjFitResponses{1}.values).^2));
    
    %% Get the time course predicitions of the CV fit CRF params to the held out data.
    %
    % First need to put the held out data into direction/contrast form.
    %
    % This conversion is possible because the IAMP packet pocket has meta data
    % that we put there to allow exactly this conversion.  That meta data
    % encapsulates the key things we need to know about the stimulus obtained
    % from the analysis parameters.
    heldOutDirectionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(heldOutTCPackets);
    
    % Get the time course predicitions from the NR fit to the CRF
    timeCoursePlot.nr = responseFromPacket('nrPred', analysisParams, nrCrfParams{1}, heldOutDirectionTimeCoursePacketPocket, 'plotColor', [0, 0, 1]);
    
    % Get the time course predicitions from the NR common amplitude fit to the CRF
    timeCoursePlot.nrAmp = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmp{1}, heldOutDirectionTimeCoursePacketPocket, 'plotColor', [.3, .2, .8]);
    
    % Get the time course predicitions from the NR common xponent fit to the CRF
    timeCoursePlot.nrExp = responseFromPacket('nrPred', analysisParams, nrCrfParamsExp{1}, heldOutDirectionTimeCoursePacketPocket, 'plotColor', [.5, .66, .9]);
    
    % Get the time course predicitions from the NR common amp and exponent fit to the CRF
    timeCoursePlot.nrAmpExp = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmpExp{1}, heldOutDirectionTimeCoursePacketPocket, 'plotColor', [0.7, .6, .1]);
    
    % Get the time course predicitions from the NR common amp, exponent and semi-saturation fit to the CRF
    timeCoursePlot.nrAmpExpSemi = responseFromPacket('nrPred', analysisParams, nrCrfParamsAmpExpSemi{1}, heldOutDirectionTimeCoursePacketPocket, 'plotColor', [0.3, .3, .1]);
    
    % Get the time course predicitions fromt the QCM params fit to the CRF
    timeCoursePlot.qcm = responseFromPacket('qcmPred', analysisParams, qcmCrfMeanParams{1}, heldOutDirectionTimeCoursePacketPocket, 'plotColor', [0, 1, 0]);
    
    % Add held out time course
    timeCoursePlot.heldOutRawTC = heldOutRawTC;
    
    % Plot the time course prediction for each run using the different fits to the crf
    plotTimeCourse(analysisParams, timeCoursePlot, heldOutBaselineShift, 2);
    
    %% Calculate RMSE to held out time course for this fold, given CV fits above
    %
    % Loop over sessions
    for kk = 1:length(timeCoursePlot.heldOutRawTC)
        
        % Calculate the RMSE for the mean model
        meanModelRMSE(kk,ii) = sqrt(mean(((cvMeanModel+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));
        
        % Calculate the RMSE for the NR preds
        nrRMSE(kk,ii) = sqrt(mean(((timeCoursePlot.nr{kk}.values+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));
        
        % Calculate the RMSE for the NR common amp preds
        nrAmpRMSE(kk,ii) = sqrt(mean(((timeCoursePlot.nrAmp{kk}.values+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));
        
        % Calculate the RMSE for the NR common exp preds
        nrExpRMSE(kk,ii) = sqrt(mean(((timeCoursePlot.nrExp{kk}.values+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));
        
        % Calculate the RMSE for the NR common amp and exp preds
        nrAmpExpRMSE(kk,ii) = sqrt(mean(((timeCoursePlot.nrAmpExp{kk}.values+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));
        
        % Calculate the RMSE for the NR common amp, exp and semi preds
        nrAmpExpSemiRMSE(kk,ii) = sqrt(mean(((timeCoursePlot.nrAmpExpSemi{kk}.values+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));
        
        % Calculate the RMSE for the qcm
        qcmRMSE(kk,ii) = sqrt(mean(((timeCoursePlot.qcm{kk}.values+heldOutBaselineShift(kk))-timeCoursePlot.heldOutRawTC{kk}.values).^2));   
    end   
end

%% Calc Mean and SEM for CV CRF RMSE
meanMmCrfRMSE = mean(mmCrfRMSE(:));
semMmCrfRMSE = std(mmCrfRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanNrCrfRMSE = mean(nrCrfRMSE(:));
semNrCrfRMSE = std(nrCrfRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanNrCrfAmpRMSE = mean(nrCrfAmpRMSE(:));
semNrCrfAmpRMSE = std(nrCrfAmpRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanNrCrfExpRMSE = mean(nrCrfExpRMSE(:));
semNrCrfExpRMSE = std(nrCrfExpRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanNrCrfAmpExpRMSE = mean(nrCrfAmpExpRMSE(:));
semNrCrfAmpExpRMSE = std(nrCrfAmpExpRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanNrCrfAmpExpSemiRMSE = mean(nrCrfAmpExpSemiRMSE(:));
semNrCrfAmpExpSemiRMSE = std(nrCrfAmpExpSemiRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanQcmCrfRMSE = mean(qcmCrfRMSE(:));
semQcmCrfRMSE = std(qcmCrfRMSE(:))./sqrt(size(heldOutRunOrder,2));

%% Plot CV RMSE results for CRF
figure;hold on
plotnames = categorical({'meanCrfMMRMSE','meanCrfNrRMSE','meanCrfNrAmpRMSE','meanCrfNrExpRMSE','meanCrfNrAmpExpRMSE','meanCrfNrAmpExpSemiRMSE','meanCrfQcmRMSE'});
barsForPlotCRF = [meanMmCrfRMSE, meanNrCrfRMSE,meanNrCrfAmpRMSE,meanNrCrfExpRMSE,meanNrCrfAmpExpRMSE,meanNrCrfAmpExpSemiRMSE,meanQcmCrfRMSE];
errorForPlotCRF = [semMmCrfRMSE, semNrCrfRMSE,semNrCrfAmpRMSE,semNrCrfExpRMSE,semNrCrfAmpExpRMSE,semNrCrfAmpExpSemiRMSE,semQcmCrfRMSE];
bar(plotnames,barsForPlotCRF)
errorbar(plotnames, barsForPlotCRF,errorForPlotCRF,'o')

%% Get the mean and SEM CV time course RMSE
meanMMRMSE = mean(meanModelRMSE(:));
semMMRMSE  = std(meanModelRMSE(:))./sqrt(size(heldOutRunOrder,2));

meanNrRMSE = mean(nrRMSE(:));
semNrRMSE = std(nrRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common amp preds
meanNrAmpRMSE = mean(nrAmpRMSE(:));
semNrAmpRMSE= std(nrAmpRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common exp preds
meanNrExpRMSE = mean(nrExpRMSE(:));
semNrExpRMSE = std(nrExpRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common amp and exp preds
meanNrAmpExpRMSE = mean(nrAmpExpRMSE(:));
semNrAmpExpRMSE = std(nrAmpExpRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the NR common amp, exp and semi preds
meanNrAmpExpSemiRMSE = mean(nrAmpExpSemiRMSE(:));
semNrAmpExpSemiRMSE = std(nrAmpExpSemiRMSE(:))./sqrt(size(heldOutRunOrder,2));

% Calculate the RMSE for the QCM
meanQcmRMSE= mean(qcmRMSE(:));
semQcmRMSE = std(qcmRMSE(:))./sqrt(size(heldOutRunOrder,2));

%% Plot CV RMSE results for time course
figure;hold on
plotnames = categorical({'meanMMRMSE','meanNrRMSE','meanNrAmpRMSE','meanNrExpRMSE','meanNrAmpExpRMSE','meanNrAmpExpSemiRMSE','meanQcmRMSE'});
barsForPlot = [meanMMRMSE, meanNrRMSE,meanNrAmpRMSE,meanNrExpRMSE,meanNrAmpExpRMSE,meanNrAmpExpSemiRMSE,meanQcmRMSE];
errorForPlot = [semMMRMSE, semNrRMSE,semNrAmpRMSE,semNrExpRMSE,semNrAmpExpRMSE,semNrAmpExpSemiRMSE,semQcmRMSE];
bar(plotnames,barsForPlot)
errorbar(plotnames, barsForPlot,errorForPlot,'o')
