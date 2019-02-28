function [fitOBJ,fitParamsCell, objFitResponses] = fitDirectionModel(analysisParams, modelType, packetPocket, varargin)
% Fits mutliple packets given a model type. 
%
% Syntax:
%   [fitOBJ,fitParamsCell, objFitResponses] = fitDirectionModel(analysisParams, modelType, packetPocket, varargin);
%
% Description:
%    This function takes in the a cell array of packets and returnd a cell
%    array of fit params, responses and the fit object. This currently
%    works for the naka-rushotn and qcm models. 
%
% Inputs:
%    analysisParams     - Struct of important information for the
%                         analysis
%    modelType          - the type of model. Either qcmFit or nrFit
%    packetPocket       - Cell array of packets to be fit
%
% Outputs:
%    fitOBJ             - The object used ofr fitting
%    fitParamsCell      - Cell array of fit params 
%    objFitResponses    - Cell array of the fit responses
%
% Optional key/value pairs:
%  'intialParams'         - Struct (default empty). Params structure
%                           containing initial parameters for search.
%                           See tfe.fitResponse for more.
%
% Optional key/value pairs for NR model fits:
%  'lockOffsetToZero'     - default false
%  'commonAmp'            - default false
%  'commonSemi'           - default false
%  'commonExp'            - default false 
%  'commonOffset'         - default true

% History:
%   MAB 09/09/18
%   MAB 01/06/19      -- changed from runIAMP_QCM to fit_IAMP and removed QCM
%   MAB, DHB 02/28/19 -- add initialParams key/value pair


%% Fit the tfeQCM to the IAMP beta weights
% allow QCM to fit the offset

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addRequired('modelType',@ischar);
p.addRequired('packetPocket',@iscell);
p.addParameter('lockOffsetToZero',false,@islogical);
p.addParameter('commonAmp',false,@islogical);
p.addParameter('commonSemi',false,@islogical);
p.addParameter('commonExp',false,@islogical);
p.addParameter('commonOffset',true,@islogical);
p.addParameter('initialParams',[],@(x)(isempty(x) | isstruct(x)));


p.parse(analysisParams,modelType,packetPocket,varargin{:});

switch modelType
    case 'qcmFit'
        clear defaultParamsInfo
        defaultParamsInfo.noOffset = false;
        fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);
        for ii = 1:length(packetPocket)
            % Fit the packet
            [fitParamsCell{ii},fVal,objFitResponses{ii}] = fitOBJ.fitResponse(packetPocket{ii},'defaultParamsInfo',defaultParamsInfo,'initialParams',p.initialParams);
            fprintf('\nQCMDirection parameters from direction fit to IAMP betas:\n');
            fitOBJ.paramPrint(fitParamsCell{ii})
        end
        
    case 'nrFit'
        
        uniqueDirections = round(analysisParams.directionCoding(1:analysisParams.theDimension,:),4);
        
        %% Fit the NRDirections to the the with the key/value pair flags
        % Create the tfeNakaRushtonDirection object
        fitOBJ = tfeNakaRushtonDirection(uniqueDirections, ...
            'lockOffsetToZero',p.Results.lockOffsetToZero,'commonAmp',p.Results.commonAmp,'commonSemi', p.Results.commonSemi, ...
            'commonExp',p.Results.commonExp,'commonOffset',p.Results.commonOffset,'initialParams',p.Results.initialParams);
        
        % loop over packets in the cell
        for ii = 1:length(packetPocket)
            % Fit the packet
            [fitParamsCell{ii},~,objFitResponses{ii}] = fitOBJ.fitResponse(packetPocket{ii},'initialParams',p.Results.initialParams);
            fprintf('\nNRDirection parameters from fit to IAMP betas with common offset:\n');
            fitOBJ.paramPrint(fitParamsCell{ii})
        end
        
end


end



