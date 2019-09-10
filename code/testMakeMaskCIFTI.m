areaNum = 1;
eccenRange = [0 90];
anglesRange = [0 180];
hemisphere = 'combined';
threshold = 0.9;

analysisParams = getSubjectParams('AP26_replication');
[~, userID] = system('whoami');
userID = strtrim(userID);

saveName = fullfile('/Users', userID, 'Desktop/lh.V1.dscalar.nii');


melaAnalysisPath = '/Users/michael/labDropbox/MELA_analysis/';
pathToBensonMasks = fullfile(melaAnalysisPath, 'mriTOMEAnalysis','flywheelOutput','benson');
pathToBensonMappingFile = fullfile(pathToBensonMasks,'indexMapping.mat');
pathToTemplateFile = fullfile(pathToBensonMasks,'template.dscalar.nii');

[ maskMatrix ] = makeMaskFromRetinoCIFTI(areaNum, eccenRange, anglesRange, hemisphere, 'saveName', saveName, 'pathToBensonMasks', pathToBensonMasks, 'pathToTemplateFile', pathToTemplateFile, 'pathToBensonMappingFile', pathToBensonMappingFile, 'threshold', threshold);

