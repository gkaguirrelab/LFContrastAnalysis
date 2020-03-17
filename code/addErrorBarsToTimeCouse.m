function [timeCoursePackets]= addErrorBarsToTimeCouse(errorVals,timeCoursePackets,varargin)
% Cuts up the error values for the full run and adds them to the time course packet. 
%
% Syntax:
%    timeCoursePackets]= addErrorBarsToTimeCouse(errorVals,timeCoursePackets,varargin)   
%
% Description:
%    Takes in a error vector and time course cell array and divides up the 
%    error values into of equal lengths and added them to each corresponding 
%    struct.
%
% Inputs:
%    errorVals                  - The vector or error values calculated for
%                                 the full experimental time course.
%    timeCoursePacket           - The time course packet cell array.
%                                 An array of structs with "timebase" and 
%                                 "values" subfeilds.  
%
% Outputs:
%    timeCoursePackets          - The updated time course cell array with
%                                 a shaddedErrorBars subfield added. 
%
% Optional key/value pairs:
%    - none for now

% MAB 03/17/20 created it

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('errorVals',@isvector);
p.addRequired('timeCoursePackets',@iscell);
p.parse(errorVals,timeCoursePackets,varargin{:});

%% Safety check
numPackets = length(timeCoursePackets);
numTimePoints = length(timeCoursePackets{1}.timebase);
if numPackets*numTimePoints ~= length(errorVals)
    error('the error values passed are not a match for the packets')
end

% Make the packets 
for ii = 1:numPackets;
    % Chop up the error values
    startPos = 1 + (ii - 1)*numTimePoints;
    stopPos = (ii)*numTimePoints;
    cutErrorVals = errorVals(startPos:stopPos);
    timeCoursePackets{ii}.shaddedErrorBars = cutErrorVals;
end