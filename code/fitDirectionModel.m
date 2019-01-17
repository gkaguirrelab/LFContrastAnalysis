function [fitOBJ,fitParamsCell, objFitResponses] = fitDirectionModel(analysisParams, modelType, packetPocket, varargin);
% Takes in the clean time series data and the analysis params and fits the IAMP model.
%
% Syntax:
%   [analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams] = fit_IAMP(analysisParams, fullCleanData);
%
% Description:
%    This function takes in the clean time series data and the analysis params
%    and fits the IMAP model. This function builds a stimulus design matirx
%    based on the analysisParams (from each run of the experiemnt) and run the
%    IAMP model on the cleaned and trial sorted data.
%
% Inputs:
%    analysisParams             - Struct of important information for the
%                                 analysis
%    fullCleanData              - The cleaned time course
%
% Outputs:
%    analysisParams             - Returns analysisParams with any updates
%    iampTimeCoursePacketPocket - Cell array of IAMP packets for each run
%    iampOBJ                    - The IAMP object
%    iampParams                 - Cell array of IAMP parameter fits for each run
%
% Optional key/value pairs:
%    none

% MAB 09/09/18
% MAB 01/06/19 -- changed from runIAMP_QCM to fit_IAMP and removed QCM

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

p.parse(analysisParams,modelType,packetPocket,varargin{:});

switch modelType
    case 'qcmFit'
        clear defaultParamsInfo
        defaultParamsInfo.noOffset = false;
        fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);
        for ii = 1:length(packetPocket)
            % Fit the packet
            [fitParamsCell{ii},fVal,objFitResponses{ii}] = fitOBJ.fitResponse(packetPocket{ii},'defaultParamsInfo',defaultParamsInfo);
            fprintf('\nQCMDirection parameters from direction fit to IAMP betas:\n');
            fitOBJ.paramPrint(fitParamsCell{ii})
        end
        
    case 'nrFit'
        
        uniqueDirections = round(analysisParams.directionCoding(1:analysisParams.theDimension,:),4);
        
        %% Fit the NRDirections to the the with the key/value pair flags
        % Create the tfeNakaRushtonDirection object
        fitOBJ = tfeNakaRushtonDirection(uniqueDirections, ...
            'lockOffsetToZero',p.Results.lockOffsetToZero,'commonAmp',p.Results.commonAmp,'commonSemi', p.Results.commonSemi, ...
            'commonExp',p.Results.commonExp,'commonOffset',p.Results.commonOffset);
        
        % loop over packets in the cell
        for ii = 1:length(packetPocket)
            % Fit the packet
            [fitParamsCell{ii},~,objFitResponses] = fitOBJ.fitResponse(packetPocket{ii});
            fprintf('\nNRDirection parameters from fit to IAMP betas with common offset:\n');
            fitOBJ.paramPrint(fitParamsCell{ii})
        end
        
end


end



