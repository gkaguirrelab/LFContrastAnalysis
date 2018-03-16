function meanSignal = extractMeanSignalFromMask(funcRuns,mask)
%extractMeanSignalFromMask -- extracts and plot the mean signal of an roi
%
% Inputs:
%   funcRuns -- BOLD timeseries nifti volume.
%   mask     -- A mask indicating the voxels of interest.
%
% Outputs:
%   meanSignal -- The mean signal across voxels of the ROI per timepoint.
%
% Key Value Pairs:
%   none
%
% Usage: 
%   meanSignal = extractMeanSignalFromMask(funcRuns,mask)

% MAB 2018 -- wrote function

% clean up ROI from the output of non linear warp
if length(unique(mask)) ~= 2
    mask(mask <0.1) = 0;
    mask(mask > 0.1 & mask <=1.0) = 1;
end

% make mask logical 
mask = logical(mask);


% extract and tkae mean of time series
for ii = 1:length(funcRuns)
    % load nifti for functional run
    nii = MRIread(funcRuns{ii});
    timeSeries = nii.vol;
    
    % index each time point (TR)
    for jj = 1:size(timeSeries,4)
        vol = timeSeries(:,:,:,jj);
        meanSignal(jj,ii) = mean(vol(mask));
    end
end

end