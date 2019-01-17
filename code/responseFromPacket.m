function [responses] =  responseFromPacket(obj, analysisParams, params, packetPocket, varargin)
% computes the response to a multiple packet/param inputs
%
% Syntax:
%    [responses] =  responseFromPacket(obj, analysisParams, params, packetPocket);
%
% Description:
%    This function takes in a fitting object, parameters, and a cell of packets
%    and returns the timecourse prediction of the model. 
%
% Inputs:
%    obj                        - Fitting objects     
%    analysisParams             - Struct of important information for the
%                                 analysis 
%    params                     - paramters of the model used for the
%                                 response calculation 
%    packetPocket               - A cell array of packets 
% Outputs:
%    responses                  - A cell array of responses  
%
% Optional key/value pairs:
%    upsampleCRF                - Upsample stimulus input for CRF. 

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('obj',@isobject);
p.addRequired('analysisParams',@isstruct);
p.addRequired('params',@isstruct);
p.addRequired('packetPocket',@iscell);
p.addParameter('upsampleCRF',false,@islogical);
p.parse(obj, analysisParams, params, packetPocket, varargin{:});


if p.Results.upsampleCRF
    contrastSpacing = linspace(max(analysisParams.contrastCoding),0,analysisParams.numSamples);
    tmpStim         = generateStimCombinations(contrastSpacing,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension);
    crfTimebase     = 1:length(tmpStim);
    [directions,contrasts] = tfeQCMStimuliToDirectionsContrasts(tmpStim,'precision',4);
    crfStimulus     = [directions;contrasts];
end

for ii = 1:length(packetPocket)
    
    if p.Results.upsampleCRF
        packetPocket{ii}.stimulus = crfStimulus;
    end
    
    responses{ii}.values   = obj.computeResponse(params,packetPocket{ii}.stimulus,packetPocket{ii}.kernel);
    
    if p.Results.upsampleCRF
        responses{ii}.timebase = crfTimebase;
    else
        responses{ii}.timebase = packetPocket{ii}.stimulus.timebase;
    end
    
    
end