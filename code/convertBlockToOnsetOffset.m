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

numPoints = 0;
if p.Results.onset
    numPoints = numPoints + 1;
end
if p.Results.midpoint
    numPoints = numPoints +1;
end
if p.Results.offset
    numPoints = numPoints +1;
end


for kk = 1:size(expParams,1)
    if expParams(kk,3) ~= baselineCondNum
        
        if p.Results.onset
            condRegMat(expParams(kk,3),expParams(kk,1),expParams(kk,4)) = 1;
        end
        
        if p.Results.midpoint
            condRegMat(expParams(kk,3),ceil((expParams(kk,1)+expParams(kk,2))./2),expParams(kk,4)) = 1;
        end
        
        if p.Results.offset
            condRegMat(expParams(kk,3),expParams(kk,2),expParams(kk,4)) = 1;
        end
        
    elseif expParams(kk,3) == baselineCondNum
        if p.Results.onset
            baselineRegVec(1,expParams(kk,1)) = 1;
        end
        
        if p.Results.midpoint
            baselineRegVec(1,ceil((expParams(kk,1)+expParams(kk,2))./2)) = 1;
        end
        
        if p.Results.offset
            baselineRegVec(1,expParams(kk,2)) = 1;
        end
        
    end
end

blockMatrix = [];
for ii = 1:length(unique(expParams(:,4)))
    blockMatrix = vertcat(blockMatrix,condRegMat(:,:,ii));
end

counter = 1;
expandedMatrix = zeros(size(blockMatrix,1)*numPoints,size(blockMatrix,2));
for jj = 1:size(blockMatrix,1)
    pos = find(blockMatrix(jj,:));
    for pp = 1:length(pos)
        expandedMatrix(counter,pos(pp)) = 1;
        counter=counter+1;
    end
end

% split up the baseline
bPos = find(baselineRegVec);
fullBaselineMatrix = zeros(numPoints,size(expandedMatrix,2));
for hh = 1:size(fullBaselineMatrix,1)
    fullBaselineMatrix(hh,bPos(hh:numPoints:end)) = 1;
end

onOffMat = vertcat(expandedMatrix,fullBaselineMatrix);

end