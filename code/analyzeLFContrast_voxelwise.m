% Initialize
%clear;

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26');

% Get the cleaned time series

switch analysisParams.preproc
    case 'fmriprep'
        [fullCleanData, analysisParams] = getTimeCourse(analysisParams);
    case 'hcp'
        [fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);
    otherwise
        error('Preprocessing method unknown')
end


% loop over voxels

for ii = 1:2%:size(fullCleanData,1)

    voxelTimeSeries= fullCleanData(ii,:,:);
    
    [qcmParams] = fitQCMtoVoxel(analysisParams,voxelTimeSeries);
    
    qcmParamsLst(ii,:) = qcmParams';
    
end