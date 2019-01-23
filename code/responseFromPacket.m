function [responses] =  responseFromPacket(obj, analysisParams, params, packetPocket, varargin)
% Compute the response to a multiple packet/param inputs
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
%
% Outputs:
%    responses                  - A cell array of responses  
%
% Optional key/value pairs:
%    none 

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('obj',@isobject);
p.addRequired('analysisParams',@isstruct);
p.addRequired('params',@isstruct);
p.addRequired('packetPocket',@iscell);
p.parse(obj, analysisParams, params, packetPocket, varargin{:});


% loop over packets 
for ii = 1:length(packetPocket)
    
    
    responses{ii}.values   = obj.computeResponse(params,packetPocket{ii}.stimulus,packetPocket{ii}.kernel);


    
    
end