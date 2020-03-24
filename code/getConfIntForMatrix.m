function [confInt] = getConfIntForMatrix(values,rowsOrColumns, varargin)
% Returns the 95% CI either per collumn or row for a matrix.
%
% Syntax:
%   [figHdl] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampsPoints);
%
% Description:
%    This function will return either the row-wise or culumn-wise
%    confidence intervals
%
% Inputs:
%    values                    - Matrix of values to want the confidence
%                                intervals for.
%    rowsOrColumns             - string. Either "row" or column" to specify
%                                the dimension of the CI.
%
% Outputs:
%    confInt                   - the row- or column-wise confidence
%                                intervals. row one = upper CI, row 2 =
%                                lower CI.
%
% Optional key/value pairs:
%    confIntRange              - Default 95% CI. Confidence interval range.

% MAB 09/09/18

% Subplot size
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('values',@ismatrix);
p.addRequired('rowsOrColumns',@isstr);
p.addParameter('confIntRange',0.95,@isnumeric);
p.parse(values,rowsOrColumns,varargin{:});

% sort the data per row or collumn
switch rowsOrColumns
    case 'row'
        sortedVals = sort(values,1);
        numElements = size(sortedVals,1);
    case 'column'
        sortedVals = sort(values,2);
        numElements = size(sortedVals,2);
end

% Find lower and upper CI index
rangeVals = (numElements.*(1-p.Results.confIntRange))./2;
upperCiIndx = numElements - rangeVals;
lowerCiIndx = rangeVals;
if lowerCiIndx <1
    lowerCiIndx = 1;
end

if round(rangeVals) == rangeVals
    switch rowsOrColumns
        case 'row'
            upperCI = sortedVals(upperCiIndx,:);
            lowerCI = sortedVals(lowerCiIndx,:);
        case 'column'
            upperCI = sortedVals(:,upperCiIndx)';
            lowerCI = sortedVals(:,lowerCiIndx)';
    end
else
    switch rowsOrColumns
        
        case 'row'
            upperCI = (sortedVals(floor(upperCiIndx),:) + sortedVals(ceil(upperCiIndx),:))./2;
            lowerCI = (sortedVals(floor(lowerCiIndx),:) + sortedVals(ceil(lowerCiIndx),:))./2;
        case 'column'
            upperCI = ((sortedVals(:,floor(upperCiIndx)) + sortedVals(:,ceil(upperCiIndx)))./2)';
            lowerCI = ((sortedVals(:,floor(lowerCiIndx)) + sortedVals(:,ceil(lowerCiIndx)))./2)';
    end
end

confInt = [upperCI;lowerCI];

end