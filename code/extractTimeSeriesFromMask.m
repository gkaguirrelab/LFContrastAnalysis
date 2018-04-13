function [voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRuns,mask,varargin)
%extractTimeseriesFromMask -- extracts the voxel values from functional runs
%                             that are defined by a binary mask 
%
% Inputs:
%   functionalRuns        - full path to BOLD timeseries nifti volume.
%   mask                  - A mask indicating the voxels of interest. This
%                           mask should be binary but if the values a lie
%                           bewteen zero and one, you can run with the
%                           'binarizeMask' true key value pair
%
% Outputs:
%   voxelTimeSeries       - A voxel by time by run matrix of the voxels from
%                           the mask.
%   voxelIndex             -The i,j,k position of voxels in the mask. This
%                           is in the same order as the the rows of the 
%                           voxelTimeSereis output. 
% Key Value Pairs:
%   binarizeMask          - Binarize the input mask is nonbinary. This is
%                           set such that any index < 0.1 get set to 0 and 
%                           >= 0.1 get set to 1. You should use this if you
%                           warp a binary mask with ANTs. 
%   threshold             - Set value of threshold
%
% Usage: 
%   [voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRuns,mask,varargin)

% MAB 03/2018 -- wrote function

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('binarizeMask',true, @islogical);
p.addParameter('threshold',0.5, @isnumeric);
p.parse(varargin{:});

% binarize mask
if (p.Results.binarizeMask)
    mask(mask <p.Results.threshold) = 0;
    mask(mask >= p.Results.threshold & mask <=1.0) = 1;
end

% assert the mask is binary
assert(all(unique(mask) == [0 1]'), 'Assertion Failed: Mask is not binary');

% make mask logical 
mask = logical(mask);

% extract and tkae mean of time series
for ii = 1:length(functionalRuns)
    % load nifti for functional run
    nii = MRIread(functionalRuns{ii});
    timeSeries = nii.vol;
    
    % index each time point (TR)
    for jj = 1:size(timeSeries,4)
        vol = timeSeries(:,:,:,jj);
        voxelTimeSeries(:,jj,ii) = vol(mask);
    end
end

%% Create maxtrix of mask voxel positions 
idx = find(mask);
[row,col,pag] = ind2sub(size(mask),idx);
voxelIndex = [row,col,pag];

end
