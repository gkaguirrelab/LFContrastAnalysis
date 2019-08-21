function [cSeries] =  censorFrames(cPoints,timeSeries)
% Returns the time points that need to be censored for a given set of motion
% estimates.
%
% Syntax:
%    [cPoints, numCPoints] = censorTimePoints(motionEstimates,varargin)
%
% Description:
%   Calculates the framewise dispalcemnt from the input motion esimates and
%   and applies a threshold to find the time points that need to be
%   censored within a given functional time series.
%
% Inputs:
%   cPoints               - Vector of frames that need to be censored.
%   timeSeries            - Any time sereies information needing to be 
%                           censored. 
%
% Outputs:
%   cPoints               - Vector of frames that need to be censored
%   percentCensored       - The percent of time frames that need to be
%                           censored
%
% Optional key/value pairs:
%   - none

% MAB 08/21/2019 -- wrote it.

%% Input Parser
p = inputParser; p.KeepUnmatched = false;
p.addRequired('cPoints', @isvector);
p.addRequired('timeSeries', @isvector);
p.parse(cPoints, timeSeries)

cSeries = timeSeries;
cSeries(cPoints) = nan;

end

