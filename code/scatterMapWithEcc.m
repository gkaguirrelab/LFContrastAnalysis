function [] = scatterMapWithEcc(analysisParams,mapOfInterest, hemi)


pathToBensonMasks = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');
pathToBensonEccFile = fullfile(pathToBensonMasks,[hemi '.benson14_angle.dscalar.nii']);
[ eccMap ] = loadCIFTI(pathToBensonEccFile);

[ mapVals ] = loadCIFTI(mapOfInterest);

eccVals = eccMap(mapVals ~= 0)

mapVals = mapVals(mapVals ~= 0)


scatter(eccVals,mapVals)


end