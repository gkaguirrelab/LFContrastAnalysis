function analyzeLFContrast_voxelwise(subjID)
% Initialize
%clear;
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjID);
sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
mapSavePath    = fullfile(sessionDir,'hcp_func','surfaceMaps');
templateFile   = '/Users/michael/labDropbox/MELA_analysis/hcpTemplates/template.dscalar.nii';

% Get the cleaned time series

[fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams);

distcomp.feature( 'LocalUseMpiexec', false )
% loop over voxels
ciftiMap = zeros(91282,7);
meanRSquaredQCMMap = zeros(91282,1);
stdRSquaredQCMMap = zeros(91282,1);
meanRSquaredIAMPMap = zeros(91282,1);
stdRSquaredIAMPMap = zeros(91282,1);
for ii = 1:size(fullCleanData,1)

    voxelTimeSeries= fullCleanData(ii,:,:);
    
    [qcmParams,meanRsquaredIAMP,stdRsquaredIAMP, meanRsquaredQCM,stdRsquaredQCM] = fitQCMtoVoxel(analysisParams,voxelTimeSeries);
    
    ciftiMap(voxelIndex(ii),:) = qcmParams';
    meanRSquaredQCMMap(voxelIndex(ii),:)  = meanRsquaredQCM;
    stdRSquaredQCMMap(voxelIndex(ii),:)   = stdRsquaredQCM;
    meanRSquaredIAMPMap(voxelIndex(ii),:) = meanRsquaredIAMP;
    stdRSquaredIAMPMap(voxelIndex(ii),:)  = stdRsquaredIAMP;

end


% write out minor axis ratio
ciftiVec = ciftiMap(:,1);
mapName        = fullfile(mapSavePath,['minorAxisMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,2);
mapName        = fullfile(mapSavePath,['angleMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,3);
mapName        = fullfile(mapSavePath,['nrAmplitudeMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,4);
mapName        = fullfile(mapSavePath,['nrSemiMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,5);
mapName        = fullfile(mapSavePath,['nrExpMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out mean R squared map
mapName        = fullfile(mapSavePath,['meanRSquaredMapQCM', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(meanRSquaredQCMMap', [], templateFile, mapName)

% write out standard devation of R squared map
mapName        = fullfile(mapSavePath,['stdRSquaredMapQCM', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(stdRSquaredQCMMap', [], templateFile, mapName)

% write out mean R squared map
mapName        = fullfile(mapSavePath,['meanRSquaredMapIAMP', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(meanRSquaredIAMPMap', [], templateFile, mapName)

% write out standard devation of R squared map
mapName        = fullfile(mapSavePath,['stdRSquaredMapIAMP', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(stdRSquaredIAMPMap', [], templateFile, mapName)
end
