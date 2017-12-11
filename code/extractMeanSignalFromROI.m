function meanSignal = extractMeanSignalFromROI(timeSeries,ROI)
%extractMeanSignalFromROI -- extracts and plot the mean signal of an roi
%
% Inputs:
%   timeSeries -- BOLD timeseries nifti volume.
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

