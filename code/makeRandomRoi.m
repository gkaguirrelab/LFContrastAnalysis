% make a random ROI not in V1-V3

pathToBensonMasks = fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');
fileName = 'combined.benson14_visualArea.dscalar.nii'
mapVals = loadCIFTI(fullfile(pathToBensonMasks,fileName));

% Initialize the map
qcmParamMap = zeros(91282,7);


saveMapName = 'V0_combined_ecc_0_to_0.dscalar.nii'