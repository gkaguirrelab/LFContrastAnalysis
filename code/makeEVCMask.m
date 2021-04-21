function [ ciftiMask ] = makeEVCMask(hemisphere, varargin)
% Make binary mask of areas map of Benson's retinotopy project.
%
% Syntax:
%  [ maskMatrix ] = makeEVCMask(hemisphere,  varargin);
%
% Description:
%  This routine makes binary retinotopy masks from Noah's project to be
%  used with CIFTI files processed through HPC's standard pipeline. 
%
% Inputs:
%  hemisphere:              - which hemisphere to be analyzed. Options include 'lh' for
% 							  left hemisphere, 'rh' for right, or 'combined' for both.
%
% Optional key-value pairs:
%  saveName					- a string which defines the full path for where to save the
%						      resulting mask. If no value is passed (the default), no mask
%							  is saved.
%  pathToBensonMasks        - a string which defines the full path to where the previously
%							  made Benson masks in CIFTI format are located.
%
% Output:
%  maskMatrix:				- a 92812 x 1 binary vector that defines the retinotopic mask.
%
%
% Example:
%{
% make a V1 mask for the left hemisphere
areaNum = 1;
eccenRange = [0 90];
anglesRange = [0 180];
hemisphere = 'lh';
threshold = 0.9;

[~, userID] = system('whoami');
userID = strtrim(userID);

saveName = fullfile('/Users', userID, 'Desktop/lh.V1.dscalar.nii');
pathToBensonMasks = fullfile('/Users', userID, 'Dropbox-Aguirre-Brainard-Lab/MELA_analysis/mriTOMEAnalysis/flywheelOutput/benson/');
pathToBensonMappingFile = fullfile('/Users', userID, 'Dropbox-Aguirre-Brainard-Lab/MELA_analysis/mriTOMEAnalysis/flywheelOutput/benson/indexMapping.mat');
pathToTemplateFile = fullfile('/Users', userID, 'Dropbox-Aguirre-Brainard-Lab/MELA_analysis/mriTOMEAnalysis/flywheelOutput/benson/template.dscalar.nii');

[ maskMatrix ] = makeEVCMask(areaNum, eccenRange, anglesRange, hemisphere, 'saveName', saveName, 'pathToBensonMasks', pathToBensonMasks, 'pathToTemplateFile', pathToTemplateFile, 'pathToBensonMappingFile', pathToBensonMappingFile, 'threshold', threshold);


% make a V1 mask for the combined right-left hemisphere
areaNum = 1;
eccenRange = [0 90];
anglesRange = [0 180];
hemisphere = 'combined';
saveName = fullfile('/Users', userID, 'Desktop/combined.V1.dscalar.nii');
pathToBensonMasks = fullfile('/Users', userID, 'Dropbox-Aguirre-Brainard-Lab/MELA_analysis/mriTOMEAnalysis/flywheelOutput/benson/');
pathToBensonMappingFile = fullfile('/Users', userID, 'Dropbox-Aguirre-Brainard-Lab/MELA_analysis/mriTOMEAnalysis/flywheelOutput/benson/indexMapping.mat');
pathToTemplateFile = fullfile('/Users', userID, 'Dropbox-Aguirre-Brainard-Lab/MELA_analysis/mriTOMEAnalysis/flywheelOutput/benson/template.dscalar.nii');

[ maskMatrix ] = makeEVCMask(areaNum, eccenRange, anglesRange, hemisphere, 'saveName', saveName, 'pathToBensonMasks', pathToBensonMasks, 'pathToTemplateFile', pathToTemplateFile, 'pathToBensonMappingFile', pathToBensonMappingFile, 'threshold', threshold);

%}

p = inputParser; p.KeepUnmatched = true;
p.addParameter('saveName', [], @ischar)
p.addParameter('pathToBensonMasks', [], @ischar)
p.addParameter('pathToBensonMappingFile', [], @ischar)
p.addParameter('pathToTemplateFile', [], @ischar)
p.parse(varargin{:});


%% Locate the template files
% describe the different templates we want to produce
mapTypes = {'angle', 'eccen', 'varea'};
hemispheres  = {'lh', 'rh'};
pathToBensonMasks = p.Results.pathToBensonMasks;

%% Restrict area

areaMask = zeros(327684,1);


if strcmp(hemisphere, 'lh') || strcmp(hemisphere, 'combined')
    rhAreaMask = zeros(163842,1);
    lhAreaMask = zeros(163842,1);
    lhAreaMap = MRIread(fullfile(p.Results.pathToBensonMasks, 'lh.benson14_varea.v4_0.mgz'));
    lhAreaMap = lhAreaMap.vol;
    
    
    
    lhAreaMask(lhAreaMap ~= 0) = 1;
    
    lhAreaMask = [lhAreaMask; rhAreaMask];
    
    areaMask = areaMask + lhAreaMask;
    
    
end
if strcmp(hemisphere, 'rh') || strcmp(hemisphere, 'combined')
    rhAreaMask = zeros(163842,1);
    lhAreaMask = zeros(163842,1);
    rhAreaMap = MRIread(fullfile(p.Results.pathToBensonMasks, 'rh.benson14_varea.v4_0.mgz'));
    rhAreaMap = rhAreaMap.vol;
    
    rhAreaMask(rhAreaMap ~= 0) = 1;
    
    rhAreaMask = [lhAreaMask; rhAreaMask];
    
    areaMask = areaMask + rhAreaMask;
    
end


%% Convert FreeSurfer mask to HCP
matrix = sparse(91282, 327684);
load(p.Results.pathToBensonMappingFile);
for ii = 1:length(ciftifsaverageix)
    matrix(ciftifsaverageix(ii), ii) = 1;
end

sumPerRow = sum(matrix,2);

nonZeroIndices = find(matrix);
for ii = 1:length(nonZeroIndices)
    [row,column] = ind2sub(size(matrix), nonZeroIndices(ii));
    matrix(row,column) = matrix(row,column)/sumPerRow(row);
end

ciftiMask = matrix * areaMask;

threshold = 0.9;

ciftiMask(ciftiMask < threshold) = 0;
ciftiMask(ciftiMask >= threshold) = 1;



% save out mask, if desired
if ~isempty(p.Results.saveName)
    makeWholeBrainMap(ciftiMask', [], fullfile(p.Results.pathToTemplateFile), p.Results.saveName)
end




end