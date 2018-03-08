function meanSignal = extractMeanSignalFromROI(funcRuns,areaMap,eccMap,area, eccThresh)
%extractMeanSignalFromROI -- extracts and plot the mean signal of an roi
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

% clean up area map from the output of non linear warp
if length(unique(areaMap)) ~= 4
    areaMap(areaMap <0.5) = 0;
    areaMap(areaMap > 0.5 & areaMap <=1.0) = 1;
    areaMap(areaMap > 1 & areaMap <=2) = 2;
    areaMap(areaMap > 2 & areaMap <=3.0) = 3;
end

% make mask
areaMap = areaMap;
areaMap(areaMap ~= area) = 0;
areaMap(areaMap == area) = 1;
areaMap = logical(areaMap);

eccMap(eccMap <=  eccThresh) = 1; 
eccMap(eccMap >  eccThresh) = 0; 
eccMap = logical(eccMap);

mask = (areaMap & eccMap) == 1;

for ii = 1:length(funcRuns)
    % load nifti for functional run
    nii = MRIread(funcRuns{ii});
    timeSeries = nii.vol;
    
    for jj = 1:size(timeSeries,4)
        vol = timeSeries(:,:,:,jj);
        meanSignal(jj,ii) = mean(vol(mask));
    end
end


% %% Debug checks
% 
% A = vol + 10000*mask;
% 
% for ii = 1:size(vol,3)
%     imagesc(A(:,:,ii)) 
%     colormap jet
%     waitforbuttonpress
% end