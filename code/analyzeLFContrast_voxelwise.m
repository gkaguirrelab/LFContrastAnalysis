% Initialize
%clear;
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('LZ23');
sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
mapSavePath    = fullfile(sessionDir,'hcp_func','surfaceMaps');
templateFile   = '/Users/michael/labDropbox/MELA_analysis/hcpTemplates/template.dscalar.nii';

% Get the cleaned time series

[fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams);

distcomp.feature( 'LocalUseMpiexec', false )
% loop over voxels
ciftiMap = zeros(91282,7);

for ii = 1:size(fullCleanData,1)

    voxelTimeSeries= fullCleanData(ii,:,:);
    
    [qcmParams,meanRsquared,stdRsquared] = fitQCMtoVoxel(analysisParams,voxelTimeSeries);
    
    ciftiMap(voxelIndex(ii),:) = qcmParams';
    meanRSquaredMap(voxelIndex(ii),:) = meanRsquared;
    stdRSquaredMap(voxelIndex(ii),:) = stdRsquared;

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
mapName        = fullfile(mapSavePath,['meanRSquaredMap', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(meanRSquaredMap', [], templateFile, mapName)

% write out standard devation of R squared map
mapName        = fullfile(mapSavePath,['stdRSquaredMap', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(stdRSquaredMap', [], templateFile, mapName)
