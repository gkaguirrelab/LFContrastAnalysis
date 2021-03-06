function [expParams] = getExpParams(dataParamFile,TR, varargin)
% getExpParams
%
% Description:
%   Takes the data file saved when OLApproach_TrialSequenceMR is run and
%   produces a matrix of block start/stop times and conditions.
%
% Inputs:
%  dataParamFile    = Full file name of the data file saved out from
%                     OLApproach_TrialSequenceMR.
%  TR               = The reption time of the scan. This is in seconds
%
% Outputs:
%   expParams       = A Nx4 matrix. Columns are (in order) 
%                     1. Stimulus block start time in TRs 
%                     2. Stimulus block end time in TRs.
%                     3. The contrast condition code. [1, .5, .25, .125, .0625 0] 
%                     4. The driection contition code.
%                     The number of rows is equal to the number of blocks 
%                     in the experiment.
%
% Optional key/value pairs:
%   stripInitialTRs = Strip the initial 2 TRs from the scan (we need this  
%                     for a data set we collected that had 2 dumt TRs in the 
%                     beginning) 
%
% Examples are provided in the source code.
%
% See also:
%

% History
%  2/07/18  mab  Created.

% Examples:
%{
    dataPath = '~/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_TrialSequenceMR/MRContrastResponseFunction/DataFiles/HERO_gka1/2017-09-19/session_1'
    fileName = 'session_1_scan3.mat';
    dataParamFile = fullfile(dataPath,fileName);
    TR = 0.800;
    expParams = getExpParams(dataParamFile,TR)
    
%}


%TR in s

p = inputParser;
p.addParameter('stripInitialTRs',true,@islogical);
p.addParameter('hrfOffset',true,@islogical);
p.parse(varargin{:});


load(dataParamFile)

%% load the stim file 

experimentStart = responseStruct.events(1).tTrialStart;

%% calculate stim duration
%this retuns a matrix that is organized in columns of 1) block start TR 2)
%block stop TR 3) contrast level 4) direction
for ii = 1: length(responseStruct.events)
    expParams(ii,1) = round((responseStruct.events(ii).tStimulusStart - experimentStart)/TR)+1;
    expParams(ii,2) = round((responseStruct.events(ii).tStimulusEnd - experimentStart)/TR);
    expParams(ii,3) = protocolParams.trialTypeOrder(1,ii);
    if size(protocolParams.trialTypeOrder,1) >=2
        expParams(ii,4) = protocolParams.trialTypeOrder(2,ii);
    end
end

if p.Results.stripInitialTRs
   expParams = expParams + [2*ones(size(expParams,1),2) zeros(size(expParams,1),1)];
end
if p.Results.hrfOffset
   expParams = expParams + [5*ones(size(expParams,1),1) zeros(size(expParams,1),2)];
end

end