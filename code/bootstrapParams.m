numBootstraps = 10;
numCores = 6;
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

parpool(numCores);

parfor ii = 1:numBootstraps
    % create the random draws with replacement
    sampleMatrix = randi([1,10],length(analysisParams.sessionFolderName),analysisParams.numAcquisitions);
    
    % random sample with replacement the iamp param fits 
    for jj  = 1:size(sampleMatrix,1)
        iampTCPacketPocketBoot(jj,:) = iampTimeCoursePacketPocket(sampleMatrix(jj,:));
        iampParamsbootstrap(jj,:) = iampParams(sampleMatrix(jj,:));
    end
    
    % get fits 
    [nrCrfParamsAmpVec(:,ii), nrCrfParamsExpVec(:,ii), nrCrfParamsAmpExpVec(:,ii), qcmCrfMeanParamsVec(:,ii)] = runDirectionModelFits(analysisParams,iampTCPacketPocketBoot,iampParamsbootstrap,iampOBJ);
end


% get the means and SEM of fits. 
% nrCrfParamsAmpMean    = 
% nrCrfParamsExpMean    =
% nrCrfParamsAmpExpMean =
% qcmCrfMeanParamsMean  = 

