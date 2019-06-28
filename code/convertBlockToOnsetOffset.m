function [onOffMat] = convertBlockToOnsetOffset(expParams,baselineCondNum,totalTime,deltaT,varargin)
% Convert a block design stimulus design matix to a maxtix that models the
% onset and offset of each block as seperate regressors.
%
% Syntax:
%   onOffMat = convertBlockToOnsetOffset(blockDesign);
%
% Description:
%    This function takes in an m x n matrix (with the rows as regressors)
%    of a block designed stimulus design matrix and converts is to to a
%    matrix that models the onset and offset of each block as seperate
%    delta function regressors.
%
% Inputs:
%    blockDesign             - An m x n matrix with the rows corresponding
%                              to binary regressors.
%
% Outputs:
%    analysisParams          - Returns analysisParams with any updates
%    baselineRowNum          - Regressors corresponding to the baseline to
%                              be modeled as a step function.
% Optional key/value pairs:
%    none


% MAB 06/04/19 wrote it.

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('expParams',@ismatrix);
p.addRequired('baselineCondNum',@isnumeric);
p.addRequired('totalTime',@isnumeric);
p.addRequired('deltaT',@isnumeric);
p.addParameter('onset',true,@islogical);
p.addParameter('midpoint',true,@islogical);
p.addParameter('offset',true,@islogical);

p.parse(expParams,baselineCondNum,totalTime,deltaT,varargin{:})
condRegMat = zeros(length(unique(expParams(:,3)))-length(baselineCondNum),totalTime/deltaT,length(unique(expParams(:,4))));
baselineRegVec = zeros(1,totalTime/deltaT);

for kk = 1:size(expParams,1)
    if expParams(kk,3) ~= baselineCondNum
        condRegMat(expParams(kk,3),expParams(kk,1):expParams(kk,2),expParams(kk,4)) = 1;
    end
    
    if p.Results.onset

    end
    
    if p.Results.midpoint
    
    end
    
    if p.Results.offset
    
    end
end
stimulusStruct.values = [];
for ii = 1:length(unique(expParams(:,4)))
    stimulusStruct.values = vertcat(stimulusStruct.values,condRegMat(:,:,ii));
end

onOffMat = vertcat(stimulusStruct.values,baselineRegVec);

end