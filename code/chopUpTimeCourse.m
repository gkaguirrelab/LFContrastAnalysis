function [thePackets] = chopUpTimeCourse(timeCoursePacket,numChops,varargin)
% Takes in a time series and chops it up into n runs.
%
% Syntax:
%    [thePackets] = chopUpTimeCourse(timeCoursePacket,numChops,varargin)   
%
% Description:
%    Takes in a time series and chops it up nto n runs of equal length.
%
% Inputs:
%    timeCoursePacket           - The time course packet to be chopped up
%                                 A struct with "timebase" and "values" 
%                                 subfeilds.  
%    numChops                   - Number of cut to be made 
%
% Outputs:
%    choppedTC                  - A cell array of the original packet 
%                                 chopped into nunChops smaller packets
%                                 with a "values" and "timebase 
% Optional key/value pairs:
%    - none for now

% MAB 12/22/19 created it

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('timeCoursePacket',@isvector);
p.addRequired('numChops',@isnumeric);
p.parse(timeCoursePacket,numChops,varargin{:});

% Cut up the time base 
numTimePoints = length(timeCoursePacket.timebase)./numChops;
if round(numTimePoints) ~= numTimePoints
    error('number of chops does not divide the time series into integer length')
end
newTimeBase = timeCoursePacket.timebase(1:numTimePoints);

% Make the packets 
for ii = 1:numChops;
    % Chop up the values
    startPos = 1 + (ii - 1)*numTimePoints;
    stopPos = (ii)*numTimePoints;
    newValues = timeCoursePacket.values(startPos:stopPos);
    
    % Make cell array of chopped packets
    thePackets{ii}.timebase = newTimeBase;
    thePackets{ii}.values = newValues;    
end