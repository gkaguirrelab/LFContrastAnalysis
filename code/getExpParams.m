function [expParams] = getExpParams(dataParamFile,TR, varargin)


%TR in s

p = inputParser;
p.addParameter('stripInitialTRs',1,@islogical);
p.parse(varargin{:});


load(dataParamFile)

%% load the stim file 

experimentStart = responseStruct.events(1).tTrialStart;

%% calculate stim duration
for ii = 1: length(responseStruct.events)
    expParams(ii,1) = round((responseStruct.events(ii).tStimulusStart - experimentStart)/TR)+1;
    expParams(ii,2) = round((responseStruct.events(ii).tStimulusEnd - experimentStart)/TR);
    expParams(ii,3) = protocolParams.trialTypeOrder(ii);
end

if stripInitialTRs
   expParams = expParams + [2*ones(size(expParams,1),2) zeros(size(expParams,1),1)];
end

end