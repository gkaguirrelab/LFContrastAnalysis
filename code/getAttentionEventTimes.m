function [attentionEventTimes eventsRegressor] = getAttentionEventTimes(block, responseStruct, varargin)
% Returns the attention event times relative to start of experiment in ms
%
% Syntax:
%   [eventsRegressor attentionEventTimes] = getAttentionEventTimes(block, responseStuct,varargin)
%
% Description:
%	Function that takes in files saved out from the OLApproach code and
%	returns a list of attention event times and an optional attention event
%	regressor. 
%
% Inputs:
%   block                 - block struct saved when running the
%                           OLApproach_TrialSequenceMR which contains
%                           attention event information 
%   responseStruct        - responcestruct saved when running the
%                           OLApproach_TrialSequenceMR which contains
%                           trial start information 
%  
% Outputs:
%   attentionEventTimes   - Vector of attention event start times relative to the
%                           start of the experiment (in ms) 
%   eventsRegressor       - A binary regressor coding the start of the
%                           attention events. The nearest TR is matchin an 
%                           attention event time is marked as 1. Returned 
%                           at the same resolution as the timebase. (only if
%                           the timebase key/value pair is provided) 
%
% Optional key/value pairs:
%   timebase              - Supply a response timebase for creation of a
%                           regressor of attention events  
%
% Examples are provided in the source code.
%

% History
%  3/30/18  mab  Created.

% Examples:
%{
    [eventsRegressor attentionEventTimes] = getAttentionEventTimes(block, responseStuct,'timebase', expTimeBase)
%}

p = inputParser; p.KeepUnmatched = false;
p.addRequired('block', @isstruct);
p.addRequired('responseStruct', @isstruct);
p.addParameter('timebase',[], @isnumeric);
p.parse(block,responseStruct, varargin{:})

% get the system start time of the experiment
expStartTime = responseStruct.tBlockStart;

% loop over the trials
count = 1;
for ii = 1:length(responseStruct.events)
    
    if  logical(block(ii).attentionTask.segmentFlag)
        
        % get the start time of the trial 
        trialStartTime = responseStruct.events(ii).tTrialStart - expStartTime;
        
        % get the time step 
        timeStep = responseStruct.timeStep;
        
        % get the attention time from start of trial  
        attentionStart =  block(ii).attentionTask.theStartBlankIndex.*timeStep + responseStruct.events(ii).trialWaitTime;
        
        % add trial start time to attention time
        attentionEventTimes(count) = trialStartTime + attentionStart;
        
        % add to count
        count = count+1;
    end
end

% convert to ms
attentionEventTimes = attentionEventTimes*1000;

% create regressor if timebase is provided
eventsRegressor = zeros(size(p.Results.timebase));
if ~isempty(p.Results.timebase)
    for jj = 1:length(attentionEventTimes)
        index = find( abs(p.Results.timebase - attentionEventTimes(jj)) == min( abs(p.Results.timebase - attentionEventTimes(jj))));
        eventsRegressor(index) = 1;
    end
end
    
end