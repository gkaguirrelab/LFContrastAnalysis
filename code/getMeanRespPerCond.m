function [] = getMeanRespPerCond(data,responseStruct,protocolParams,varargin)



p = inputParser;
p.addParameter('stripInitialTRs',1,@islogical);
p.addParameter('dimensions',3,@isnumeric);
p.parse(varargin{:});

protocolParams.trialTypeOrder

%% load the stim file 


%% calculate stim duration
for ii = 1: length(responseStruct.events)
    blockTimes(ii) = (responseStruct.events(ii).tStimulusEnd - responseStruct.events(ii).tStimulusStart);
end
%% convert to tr index
blockTRs = round(blockTimes./(func.tr/1000));

%% bin 


end