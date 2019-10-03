% Initialize
%clear;
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26');
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
    
    qcmParams = fitQCMtoVoxel(analysisParams,voxelTimeSeries);
    
    ciftiMap(voxelIndex(ii),:) = qcmParams';


end


ciftiMap(voxelIndex,:) = qcmParams;

% write out minor axis ratio
ciftiVec = ciftiMap(:,1);
mapName        = fullfile(mapSavePath,'minorAxisMap.dscalar.nii');
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,2);
mapName        = fullfile(mapSavePath,'angleMap.dscalar.nii');
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,3);
mapName        = fullfile(mapSavePath,'nrAmplitudeMap.dscalar.nii');
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,4);
mapName        = fullfile(mapSavePath,'nrSemiMap.dscalar.nii');
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = ciftiMap(:,5);
mapName        = fullfile(mapSavePath,'nrExpMap.dscalar.nii');
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)