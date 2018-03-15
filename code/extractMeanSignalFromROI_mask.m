function meanSignal = extractMeanSignalFromROI_mask(funcRuns,mask)
%extractMeanSignalFromROI_mask -- extracts and plot the mean signal of an roi
%
% Inputs:
%   funcRuns -- BOLD timeseries nifti volume.
%   ROI        -- The nifti file that tage the voxels of interest.
%
% Outputs:
%   meanSignal -- The mean signal across voxels of the ROI per timepoint.
%
% Key Value Pairs:
%   verbose    -- Plot the time course.
%
% Usage:
%   meanSignal = extractMeanSignalFromROI(timeSeries,ROI,'verbose',false)

% MAB 2017 -- started function

% clean up ROI from the output of non linear warp
if length(unique(mask)) ~= 2
    mask(mask <0.1) = 0;
    mask(mask > 0.1 & mask <=1.0) = 1;
end

% make mask

mask = logical(mask);

for ii = 1:length(funcRuns)
    % load nifti for functional run
    nii = MRIread(funcRuns{ii});
    timeSeries = nii.vol;
    
    for jj = 1:size(timeSeries,4)
        vol = timeSeries(:,:,:,jj);
        meanSignal(jj,ii) = mean(vol(mask));
    end
end

end