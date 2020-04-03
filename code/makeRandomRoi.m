% make a random ROI not in V1-V3

% number of veticies in the ROI
numVertex = 860; 

% load benson atlas
pathToBensonMasks = fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');
fileName = 'combined.benson14_visualArea.dscalar.nii'
mapVals = loadCIFTI(fullfile(pathToBensonMasks,fileName));

% Non retino areas 
nonRetinoIndx = find(mapVals == 0);
numPts = length(nonRetinoIndx);
pts = sort(randperm(numPts,numVertex));
roiPts = nonRetinoIndx(pts);

% Initialize the map
qcmParamMap = zeros(91282,1);
qcmParamMap(roiPts) = 1;

mapSavePath = fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'),'LFContrastAnalysis','MNI_ROIs');
saveMapName = fullfile(mapSavePath,'V0_combined_ecc_0_to_0.dscalar.nii');

%save it 
dropBoxPath     = fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'),'LFContrastAnalysis');
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');
makeWholeBrainMap(qcmParamMap', [], templateFile, saveMapName);