function [crfStimulus] = upsampleCRF(analysisParams,varargin)
% Provides a higher resolution contrast/direcetions base for CRF predictions
% 
% Syntax:
%   [crfStimulus] = upsampleCRF(analysisParams)
%             
% Description:
%   This function takes in a fitting object, parameters, and a cell of packets
%   and returns the timecourse prediction of the model. 
%
% Inputs:    
%   analysisParams      - Struct of important information for the
%                         analysis. Relevant fields for this are:
%                           * contrastCoding 
%                           * directionCoding 
%                           * maxContrastPerDirection 
%                           * theDimention 
%                           * numSamples - upsample resolution 
% Outputs:
%   crfStimulus         - Upsampled contrast/directions stimuli. 
%
% Optional key/value pairs:
%   none

% History:
%   01/22/2019 MAB Wrote it. 

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.parse(analysisParams, varargin{:});

% upsample the contrast spacing
contrastSpacing = linspace(max(analysisParams.contrastCoding),0,analysisParams.numSamples);

% all combinations of directions and contrast spacing (scaled by the max
% contrast per directions 
tmpStim         = generateStimCombinations(contrastSpacing,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension);

% convert tmpStim to directions and contrasts 
[directions,contrasts] = tfeQCMStimuliToDirectionsContrasts(tmpStim,'precision',4);

% match format expected by tfe
crfStimulus.values     = [directions;contrasts];
crfStimulus.timebase   = 1:size(crfStimulus.values,2);

end