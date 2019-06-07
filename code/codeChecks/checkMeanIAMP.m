% Check mean IAMP fits


%% Get the parameters
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26_replication');

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

%% Get the cleaned time series
[fullCleanData, analysisParams] = getTimeCourse(analysisParams);




%% Fit the IAMP to each run
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampRegFitResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData,'modelOnOff', false, 'plotColor', [.8,.4,.1]);


% Concatenate the runs and the stim design matrix 
%% Fit the IAMP to each run
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData,'modelOnOff', true, 'plotColor', [.3,.7,.4]);



% 
% Get the offset
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end


%% Plotting
iampPlot.timeCourse = rawTC;
iampPlot.iampFit    = iampRespiampRegFitResponsesonses;

tcHndl = plotTimeCourse(analysisParams, iampPlot, concatBaselineShift, analysisParams.numSessions*analysisParams.numAcquisitions);

