function [onOffMat] = convertBlockToOnsetOffset(blockDesign,baselineRowNum)
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

for ii = 1:size(blockDesign,1)
    
    if ii ~= baselineRowNum % exclude baseline condition
        
        % Find where onsets and offsets occur
        diffVec = diff(blockDesign(ii,:));
        
        if blockDesign(ii,1) == 1
            diffVec = [1 diffVec];
        else
            diffVec = [0 diffVec];
        end
        
        onVec = double(diffVec == 1);
        offVec = double(diffVec == -1);
        
        onOffMat((2*ii-1),:) = onVec;
        onOffMat((2*ii),:) = offVec;
        
    end
    
end

% add last trial offset
% find the empty regressor and and a one to the last postition

onOffMatSum = sum(onOffMat,2);

onOffMat(find(onOffMatSum == 0),end) = 1;

% Add the baseline as seperate blocks

% pull out  the baseline regressor
baselineCond = blockDesign(baselineRowNum,:);

% Get onset and offset postitions
baseOnOff = diff(baselineCond);
if baselineCond(1) == 1
    baseOnOff = [1 baseOnOff];
else
    baseOnOff = [0 baseOnOff];
end

if baselineCond(end) == 1
    baseOnOff(end) =-1;
end

% find the start of each baseline block
startPos = find(baseOnOff == 1);
stopPos  = find(baseOnOff == -1);

% Check for back to back baseline blocks

durations = stopPos - startPos;

% find the double block if present account for some blocks showing at 1 less
% TR
doubleBlockIndx = find(durations./min(durations) >= 1.9);

if ~isempty(doubleBlockIndx)
    for jj = 1:length(doubleBlockIndx)
        splitBlockTime = (stopPos(doubleBlockIndx(jj)) - startPos(doubleBlockIndx(jj)))./2;
        startPos = round([startPos(1:doubleBlockIndx(jj)), startPos(doubleBlockIndx(jj))+splitBlockTime, startPos(doubleBlockIndx(jj)+1:end)]);
        stopPos = round([stopPos(1:doubleBlockIndx(jj)-1) , stopPos(doubleBlockIndx(jj))-splitBlockTime, stopPos(doubleBlockIndx(jj):end)]);
    end
end

baseOnOffMat = zeros(length(startPos),size(blockDesign,2));

for kk = 1:length(startPos)
    baseOnOffMat(kk,startPos(kk):stopPos(kk)) = 1;
end

onOffMat = [onOffMat;baseOnOffMat];


