function [figHndl] = scatterMapWithEcc(subjId, mapOfInterest, varargin)
% Creates a scatter plot of a parameter of interest and eccentricity
%
% Syntax:
%   [figHndl] = scatterMapWithEcc(analysisParams,mapOfInterest, hemi)
%
% Description:
%    This function takes in a cell array of packets and returns one packet
%    with a the fields being a concatenation of all the input packets.
%
% Inputs:
%    subjId                     - String with subject ID
%    mapOfInterest              - String with map name. Either:
%                                 'minorAxis'
%                                 'angle'
%                                 'amplitude'
%                                 'semi'
%                                 'exponent'
%                                 'rSqaured'
% Outputs:
%    figHndl                    - Figure Handle
%
% Optional key/value pairs:
%    saveFigs                   - Logical flag to save the figure

% MAB 03/10/20 Wrote it.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('subjId',@ischar);
p.addRequired('mapOfInterest',@ischar);
p.addParameter('saveFigs',true,@islogical)
p.addParameter('dotColor',[1 0.5 0.5],@isvector)
p.parse(subjId,mapOfInterest,varargin{:});

% load subject params
analysisParams = getSubjectParams(subjId);

% This is sent the getSubjectParams function
hemi  = analysisParams.hemisphere;

%% LOAD THE PARAMTER MAP
% Path the subject map data
dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
mapSavePath    = fullfile(dropBoxPath,'surfaceMaps',analysisParams.expSubjID,'V1');

switch mapOfInterest
    case 'minorAxis'
        subjMapName = fullfile(mapSavePath,['minorAxisMap_', analysisParams.sessionNickname '.dscalar.nii']);
    case 'angle'
        subjMapName = fullfile(mapSavePath,['angleMap_', analysisParams.sessionNickname '.dscalar.nii']);
    case 'amplitude'
        subjMapName = fullfile(mapSavePath,['nlAmpMap_', analysisParams.sessionNickname '.dscalar.nii']);
    case 'semi'
        subjMapName = fullfile(mapSavePath,['nlSemiMap_', analysisParams.sessionNickname '.dscalar.nii']);
    case 'exponent'
        subjMapName = fullfile(mapSavePath,['nlExpMap_', analysisParams.sessionNickname '.dscalar.nii']);
    case 'rSqaured'
        subjMapName = fullfile(mapSavePath,['rSquaredMapQcm_', analysisParams.sessionNickname '.dscalar.nii']);
end

% Load the paramter map
mapVals = loadCIFTI(subjMapName);

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
voxelIndex = find(maskMatrix);

%% APPLY THE MASK
eccSatterPoints    = eccMap(voxelIndex);
paramScatterPoints = mapVals(voxelIndex);


%% PLOT IT
figHndl = figure;
hold on;
% Set the figure's size in inches
figureSizeInches = [12 8];
figHndl.Units = 'inches';

markerSize = 9;
markerAreaPtsSquared = markerSize^2;

if strcmp(mapOfInterest, 'minorAxis') | strcmp(mapOfInterest, 'angle')
    
    % Load the R^2 map
    r2MapName = fullfile(mapSavePath,['rSquaredMapQcm_', analysisParams.sessionNickname '.dscalar.nii']);
    r2Vals = loadCIFTI(r2MapName);
    r2SatterPoints    = r2Vals(voxelIndex);
    
    % Convert R^2 to color scale
    r2AlphaVals = r2SatterPoints./max(r2SatterPoints);
    
    mdlr = fitlm(eccSatterPoints,paramScatterPoints,'RobustOpts','on');
    regLineParams = mdlr.Coefficients.Variables;
    regLine = @(x) regLineParams(2,1).*x + regLineParams(1,1);
    
    for ii = 1:length(eccSatterPoints)
        scttrPltHndl= scatter(eccSatterPoints(ii),paramScatterPoints(ii), markerAreaPtsSquared, 'o', ...
            'LineWidth', 1.0, 'MarkerFaceColor', p.Results.dotColor, 'MarkerEdgeColor', 1-((r2AlphaVals(ii).*[.9 .9 .9])+0.1));
        set(scttrPltHndl, 'MarkerFaceAlpha', r2AlphaVals(ii));
    end
    xPts = [min(eccSatterPoints), max(eccSatterPoints)];
    yPts = regLine(xPts);
    l1 = line(xPts,yPts,'Color',p.Results.dotColor,'LineWidth',2);
    
else
    
    scttrPltHndl= scatter(eccSatterPoints,paramScatterPoints, markerAreaPtsSquared, 'o', ...
        'LineWidth', 1.0, 'MarkerFaceColor', p.Results.dotColor, 'MarkerEdgeColor', p.Results.dotColor);
    set(scttrPltHndl, 'MarkerFaceAlpha', 0.6);
end

% add text
modelTxtTheta = sprintf('slope = %s offset = %s',...
                num2str(regLineParams(2,1),3), num2str(regLineParams(1,1),3));
theTextHandle = text(gca, 1,.9 , modelTxtTheta, 'Interpreter', 'latex');
set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);

xlabel('Eccentricity (Degrees)');
set(figHndl, 'Renderer', 'Painters');
switch mapOfInterest
    case 'minorAxis'
        yString = 'Minor Axis Ratio';
        ylim([0 1]);
    case 'angle'
        yString = 'Ellipse Angle (Degrees)';
        ylim([-100 100]);
    case 'amplitude'
        yString = 'Non-Linearity Amplitude';
    case 'semi'
        yString = 'Non-Linearity Semisaturation';
    case 'exponent'
        yString = 'Non-Linearity Exponent' ;
    case 'rSqaured'
        yString = 'R Sqaured';
        ylim([0 1]);
end

ylabel(yString);

title([yString ' By Eccentricity']);

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

if p.Results.saveFigs
    set(figHndl, 'PaperSize',figureSizeInches);
    set(figHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
    % Full file name
    figName =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
        [analysisParams.expSubjID,'_scatter_' mapOfInterest '_' analysisParams.sessionNickname '_hcp.pdf']);
    % Save it
    print(figHndl, figName, '-dpdf', '-r300');
end
end