function [theModelPreds] =  responseFromPacket(predictionsType, analysisParams, fitParams, packetPocket, varargin)
% Compute the response to a multiple packet/param inputs
%
% Syntax:
%    [theModelPreds] =  responseFromPacket(predictionsType, analysisParams, fitParams, packetPocket, varargin)%
% Description:
%    This function takes in a fitting object, parameters, and a cell of packets
%    and returns the timecourse prediction of the model.
%
% Inputs:
%    predictionsType            - Model used for predicition (either nrPred
%                                 or qcmPred)
%    analysisParams             - Struct of important information for the
%                                 analysis
%    fitParams                     - paramters of the model used for the
%                                 response calculation
%    packetPocket               - A cell array of packets
%
% Outputs:
%    theModelPreds              - A cell array of responses
%
% Optional key/value pairs:
%    plotColor                  - rgb vector setting the color of the line

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('predictionsType',@ischar);
p.addRequired('analysisParams',@isstruct);
p.addRequired('params',@isstruct);
p.addRequired('packetPocket',@iscell);
p.addParameter('plotColor',[1 , 1, 1],@isvector);
p.parse(predictionsType, analysisParams, fitParams, packetPocket, varargin{:});


% loop over packets
for ii = 1:size(packetPocket,1)
    
    switch predictionsType
        case 'nrPred'
            srt = analysisParams.numDirPerSession*(ii - 1) + 1;
            stp = analysisParams.numDirPerSession*ii;
            directions = analysisParams.directionCoding(1:analysisParams.theDimension, srt:stp);
            fitOBJ = tfeNakaRushtonDirection(directions);
            
            params = fitParams(srt:stp);
        case 'qcmPred'
            fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);
            params = fitParams;
            
        case 'IAMP'
            fitOBJ = tfeIAMP('verbosity','none');
            if ii <=10
                params = fitParams.sessionOne;
            else
                params = fitParams.sessionTwo;
            end
        otherwise
            error('Model not known');
    end
    
    for jj = 1:size(packetPocket,2)
        
        theModelPreds{ii,jj} = fitOBJ.computeResponse(params,packetPocket{ii,jj}.stimulus,packetPocket{ii,jj}.kernel);
        theModelPreds{ii,jj}.plotColor   = p.Results.plotColor;
        
    end
end