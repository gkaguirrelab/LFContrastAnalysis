function [regLineParams] =  analyzeLFContrast_GLM_HighLow(subjId)

% High/Log Contrast change with ecc.

display(['STARTING - Main Analysis: ',subjId])
% Load the subject relevant info
analysisParams = getSubjectParams(subjId);

analysisParams.preproc  = 'hcp';
analysisParams.saveMaps = false;
analysisParams.showFigs = false;

hemi = 'combined';
sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
mapSavePath    = fullfile(dropBoxPath,'surfaceMaps',analysisParams.expSubjID);
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');

%% Set up stuff
% Initialize the maps
lMinusM_highContrast = zeros(91282,1);
lMinusM_lowContrast  = zeros(91282,1);
lPlusM_highContrast  = zeros(91282,1);
lPlusM_lowContrast   = zeros(91282,1);
% bandpass the signal
analysisParams.highpass = false;

%% Load the relevant data (SDM, HRF, TC)

%set the HRF
[analysisParams] = loadHRF(analysisParams);

if analysisParams.highpass
    analysisParams.HRF.values = highpass(analysisParams.HRF.values ,5/288,1/.8);
end

%% Get the cleaned time series
[fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams);

% reshape the data to voxels x time point(all 20 runs)
timeCourses = [];
for ii = 1:size(fullCleanData,3)
    timeCourses = [timeCourses,fullCleanData(:,:,ii)];
end

% Get a packet for each run (1-20)
[analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData,'highpass',analysisParams.highpass);

% Concatenate Packet
[analysisParams, theFullPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);

% Create a new packet
% Code a new stimulus design matrix for L-M and L+M high/low contrast
% new matrix is [L-M(High);L-M(Low);L+M(High);L-M(Low)];

numTimePoints = analysisParams.expLengthTR*analysisParams.numAcquisitions;

thePacket.stimulus.values = [sum(theFullPacket.stimulus.values(1:2,1:numTimePoints));sum(theFullPacket.stimulus.values(3:4,1:numTimePoints));...
    sum(theFullPacket.stimulus.values(6:7,1:numTimePoints));sum(theFullPacket.stimulus.values(8:9,1:numTimePoints));...
    theFullPacket.stimulus.values(41,1:numTimePoints)];
thePacket.stimulus.timebase =  theFullPacket.stimulus.timebase(1:numTimePoints);

thePacket.response.values =theFullPacket.response.values(1:numTimePoints);
thePacket.response.timebase=theFullPacket.response.timebase(1:numTimePoints);

thePacket.kernel.values   = [];
thePacket.kernel.timebase = theFullPacket.kernel.timebase(1:numTimePoints);

thePacket.metaData = theFullPacket.metaData;


%% FIT THE TIME COURSE
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');
% fit the IAMP model
defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);


%% loop over voxels
for ii = 1:size(fullCleanData,1)
    % add the voxel time course to the packet
    thePacket.response.values = timeCourses(ii,1:numTimePoints);
    
    % fit the IAMP model
    defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);
    [iampParams,fVal,iampResponses] = iampOBJ.fitResponse(thePacket,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    % subtract the baseline
    theParams = iampParams.paramMainMatrix - iampParams.paramMainMatrix(end);
    
    % Put values back in map
    lMinusM_highContrast(voxelIndex(ii),:) = theParams(1);
    lMinusM_lowContrast(voxelIndex(ii),:)  = theParams(2);
    lPlusM_highContrast(voxelIndex(ii),:)  = theParams(3);
    lPlusM_lowContrast(voxelIndex(ii),:)   = theParams(4);
end

if analysisParams.saveMaps
    % write out minor axis ratio
    mapName = fullfile(mapSavePath,['lMinusM_highContrast_', analysisParams.sessionNickname '.dscalar.nii']);
    makeWholeBrainMap(lMinusM_highContrast', [], templateFile, mapName)
    
    % write out minor axis ratio
    mapName = fullfile(mapSavePath,['lMinusM_lowContrast_', analysisParams.sessionNickname '.dscalar.nii']);
    makeWholeBrainMap(lMinusM_lowContrast', [], templateFile, mapName)
    
    % write out minor axis ratio
    mapName = fullfile(mapSavePath,['lPlusM_highContrast_', analysisParams.sessionNickname '.dscalar.nii']);
    makeWholeBrainMap(lPlusM_highContrast', [], templateFile, mapName)
    
    % write out minor axis ratio
    mapName = fullfile(mapSavePath,['lPlusM_lowContrast_', analysisParams.sessionNickname '.dscalar.nii']);
    makeWholeBrainMap(lPlusM_lowContrast', [], templateFile, mapName)
end

%% LOAD THE RETINO MAPS
% Path the benson atlas maps

pathToBensonMasks = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');

if strcmp(hemi, 'lh')
    pathToLHLhEccFile = fullfile(pathToBensonMasks,'lh.benson14_eccen.dscalar.nii');
    eccMap = loadCIFTI(pathToLHLhEccFile);
    
elseif strcmp(hemi, 'rh')
    pathToRhEccFile = fullfile(pathToBensonMasks, 'rh.benson14_eccen.dscalar.nii');
    eccMap = loadCIFTI(pathToRhEccFile);
    
elseif strcmp(hemi, 'combined')
    % Load lh and rh maps
    pathToLHLhEccFile = fullfile(pathToBensonMasks,'lh.benson14_eccen.dscalar.nii');
    lhEccMap = loadCIFTI(pathToLHLhEccFile);
    pathToRhEccFile = fullfile(pathToBensonMasks,'rh.benson14_eccen.dscalar.nii');
    rhEccMap = loadCIFTI(pathToRhEccFile);
    
    % convert nan to 0 and combine maps
    lhEccMap(find(isnan(lhEccMap))) = 0;
    rhEccMap(find(isnan(rhEccMap))) = 0;
    eccMap = lhEccMap +rhEccMap;
else
    error('Map type not found');
end

%% LOAD THE MASK USED IN THE ANALYSIS
maskRoiPath  = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis','MNI_ROIs');
maskName     = ['V', num2str(analysisParams.areaNum), '_', analysisParams.hemisphere, '_ecc_', num2str(analysisParams.eccenRange(1)), '_to_', num2str(analysisParams.eccenRange(2)),'.dscalar.nii'];
maskFullFile = fullfile(maskRoiPath,maskName);
maskMatrix = loadCIFTI(maskFullFile);
maskIndex = find(maskMatrix);

eccSatterPoints    = eccMap(maskIndex);

%% PLOT IT

% plotting color
lMinusMColorDot  = [215,59,62]./255;
lMinusMColorLine = [244,194,194]./255;
lPlusMColorDot   = [42,82,190]./255;
lPlusMColorLine  = [204,204,255]./255;

% Set the figure's size in inches
figureSizeInches = [20 12];
figHndl.Units = 'inches';

markerSize = 7;
markerAreaPtsSquared = markerSize^2;


mdlr = fitlm(eccSatterPoints,lMinusM_highContrast(maskIndex),'RobustOpts','off');
regLineParams.high.lMinusM = mdlr.Coefficients.Variables;
regLine = @(x) regLineParams.high.lMinusM(2,1).*x + regLineParams.high.lMinusM(1,1);
xPts = [min(eccSatterPoints), max(eccSatterPoints)];
yPts1 = regLine(xPts);


mdlr = fitlm(eccSatterPoints,lPlusM_highContrast(maskIndex),'RobustOpts','off');
regLineParams.high.lPlusM = mdlr.Coefficients.Variables;
regLine = @(x) regLineParams.high.lPlusM(2,1).*x + regLineParams.high.lPlusM(1,1);
yPts2 = regLine(xPts);



if analysisParams.showFigs
    
    subplot(1,2,1);
    hold on;
    axis square
    scttrPlt1= scatter(eccSatterPoints,lMinusM_highContrast(maskIndex), markerAreaPtsSquared, 'o', ...
        'LineWidth', 1.0, 'MarkerFaceColor', lMinusMColorDot, 'MarkerEdgeColor', lMinusMColorDot);
    set(scttrPlt1, 'MarkerFaceAlpha', 0.6);
    
    scttrPlt2= scatter(eccSatterPoints,lPlusM_highContrast(maskIndex), markerAreaPtsSquared, 'o', ...
        'LineWidth', 1.0, 'MarkerFaceColor', lPlusMColorDot, 'MarkerEdgeColor', lPlusMColorDot);
    set(scttrPlt2, 'MarkerFaceAlpha', 0.6);
    l1 = line(xPts,yPts1,'Color',lMinusMColorLine,'LineWidth',3);
    l2 = line(xPts,yPts2,'Color',lPlusMColorLine,'LineWidth',3);
    
    ylim([-0.5,1.0])
    xlabel('Eccentricity')
    ylabel('GLM Beta Weight')
    title('High Contrast Condition');
    set(gca, ...
        'XColor', [0.2 0.2 0.2], ...
        'YColor', [0.2 0.2 0.2], ...
        'FontName', 'Helvetica', ...
        'FontSize', 14, ...
        'FontWeight', 'normal', ...
        'TickLength',[0.01 0.01], ...
        'TickDir', 'out', ...
        'LineWidth', 0.7, ...
        'Box', 'off');
end

mdlr = fitlm(eccSatterPoints,lMinusM_lowContrast(maskIndex),'RobustOpts','off');
regLineParams.low.lMinusM = mdlr.Coefficients.Variables;
regLine = @(x) regLineParams.low.lMinusM(2,1).*x + regLineParams.low.lMinusM(1,1);
yPts3 = regLine(xPts);


mdlr = fitlm(eccSatterPoints,lPlusM_lowContrast(maskIndex),'RobustOpts','off');
regLineParams.low.lPlusM = mdlr.Coefficients.Variables;
regLine = @(x) regLineParams.low.lPlusM(2,1).*x + regLineParams.low.lPlusM(1,1);
yPts4 = regLine(xPts);


if analysisParams.showFigs
    subplot(1,2,2);
    axis square
    hold on;
    scttrPlt3= scatter(eccSatterPoints,lMinusM_lowContrast(maskIndex), markerAreaPtsSquared, 'o', ...
        'LineWidth', 1.0, 'MarkerFaceColor', lMinusMColorDot, 'MarkerEdgeColor', lMinusMColorDot);
    set(scttrPlt3, 'MarkerFaceAlpha', 0.6);
    
    scttrPlt4= scatter(eccSatterPoints,lPlusM_lowContrast(maskIndex), markerAreaPtsSquared, 'o', ...
        'LineWidth', 1.0, 'MarkerFaceColor', lPlusMColorDot, 'MarkerEdgeColor', lPlusMColorDot);
    set(scttrPlt4, 'MarkerFaceAlpha', 0.6);
    
    l3 = line(xPts,yPts3,'Color',lMinusMColorLine,'LineWidth',3);
    l4 = line(xPts,yPts4,'Color',lPlusMColorLine,'LineWidth',3);
    
    ylim([-0.5,1.0])
    xlabel('Eccentricity')
    ylabel('GLM Beta Weight')
    title('Low Contrast Condition');
    set(gca, ...
        'XColor', [0.2 0.2 0.2], ...
        'YColor', [0.2 0.2 0.2], ...
        'FontName', 'Helvetica', ...
        'FontSize', 14, ...
        'FontWeight', 'normal', ...
        'TickLength',[0.01 0.01], ...
        'TickDir', 'out', ...
        'LineWidth', 0.7, ...
        'Box', 'off');
    legend([scttrPlt3,scttrPlt4],{'L-M', 'L+M'})
    
    set(figHndl, 'PaperSize',figureSizeInches);
    set(figHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
    % Full file name
    figName =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_scatter_High_Low_Contrast_GLM' analysisParams.sessionNickname '_hcp.pdf']);
    % Save it
    print(figHndl, figName, '-dpdf', '-r300');
end

display(['COMPLETED: ',subjId])


end



