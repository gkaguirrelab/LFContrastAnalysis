function [choppedTC] = chopUpTimeCourse(timecourse,numChops,varargin)
% Takes in a time series and chops it up into n runs.
%
% Syntax:
%    [choppedTC] = chopUpTimeCourse(timecourse,numChops,varargin)   
%
% Description:
%    Takes in a time series and chops it up nto n runs of equal length.
%
% Inputs:
%    timecourse                 - The time course to be chopped up
%    numChops                   - Number of cut to be made 
%
% Outputs:
%    choppedTC                  - The the original TC chopped into n run
%                                 determied by numChops
%
% Optional key/value pairs:
%    - none for now

% MAB 12/22/19 created it

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('timecourse',@isvector);
p.addRequired('numChops',@isnumeric);
p.parse(timecourse,numChops,varargin{:});