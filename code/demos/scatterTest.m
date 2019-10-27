
analysisParams = getSubjectParams('AP26');

hemi          = 'lh' 
mapOfInterest = 'minorAxisMap' 


pathToMask = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'), 'LFContrastAnalysis','surfaceMaps',analysisParams.expSubjID);
minorAxisFullFile = fullfile(pathToMask,[mapOfInterest '_' analysisParams.sessionNickname '.dscalar.nii']);

rSquaredFullFile = fullfile(pathToMask,['meanRSquaredMapQCM' analysisParams.sessionNickname '.dscalar.nii']);



pathToBensonMasks = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');
pathToBensonEccFile = fullfile(pathToBensonMasks,[hemi '.benson14_angle.dscalar.nii']);
[ eccMap ] = loadCIFTI(pathToBensonEccFile);

[ minorAxisVals ] = loadCIFTI(minorAxisFullFile);
[ rSquaredVals ] = loadCIFTI(rSquaredFullFile);


eccVals = eccMap(minorAxisVals ~= 0);
rSquaredVals = rSquaredVals(minorAxisVals ~= 0);
minorAxisVals = minorAxisVals(minorAxisVals ~= 0);

figure
scatter3(eccVals,minorAxisVals,rSquaredVals)
xlabel('eccentricity')
ylabel('minor axis vals')
zlabel('R^2 values')
axis square




%scatterMapWithEcc(analysisParams,mapFullFile, hemi)