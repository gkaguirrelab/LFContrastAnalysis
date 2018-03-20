function [avgPerCond, conditions] = sortDataByConditions(meanSignal,expParams)
% sortDataByConditions
%
% Description:
%   Takes the timecourse data and protocol file saved when OLApproach_TrialSequenceMR
%   is run and plots the power levels, the raw average timecourse and the percent signal 
%   change relative to the mean.
%
% Inputs:
%   meanSignal      = The TR (mean of all voxels in ROI) by run matrix from 
%                     extractMeanSignalFromMask.m 
%   expParams       = The information about the trial order of the exp
%                     from loading the param file save from each run of 
%                     OLApproach_TrialSequenceMR
%
% Outputs:
%   avgPerCond      = A contition by run matrix where each cell it the mean
%                     of all time points of a particular condition. 
%   conditions      = The order of conditions (rows) of the avgPerCond martix. 
%
% Optional key/value pairs:
%   none
%
% Example: 
%   [avgPerCond, conditions] = sortDataByConditions(meanSignal,expParams)

% History
%  3/18  mab  Created.

for ii = 1:length(expParams)
    blockAvg(ii) = mean(meanSignal(expParams(ii,1):expParams(ii,2)));
end

conditions = unique(expParams(:,3));

for jj = 1:length(conditions)
    avgPerCond(jj) = mean(blockAvg(expParams(:,3) == conditions(jj))); 
end

end