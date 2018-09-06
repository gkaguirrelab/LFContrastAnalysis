function stimCombinations = generateStimCombinations(contrastSpacing,directionCoding,maxContrastPerDir,theDimension)
% Returns a matrix that contains the stimulus combinations based on the inputs. 
% Syntax:
%   stimOrder = generateStimOrder(contrastCoding,directionCoding,maxContrastPerDir)
%
% Description:
%    Returns a matrix that is the total combinations of each contrast level 
%    with each direction scaled to the maximum contrast for that direction.  
%    Each column of the output matrix is single stim.  
%
% Inputs:
%    contrastCoding    - the contrast level spacing to be applied to all
%                        directions 
%    directionCoding   - A 3xn (number of directions) matrix where the column are 
%                        specifying the weighting on each cone. The convetion is row 
%                        1 is L, 2 is M , and 3 is S. ex. [1;1;0] is L+M.
%    maxContrastPerDir - Maximum contrast in each direction. 
%    theDimension      - The dimensionality of stimulus. ex. 2 for L and M
%                        directed stimulus  and 3 for L,M, and S. 
%
% Outputs:
%    stimCombinations  - The contrast level on each L,M,S for each
%    diecrtion and contrast level (scaled to the max contrast per
%    direction).
%
% Optional key/value pairs:
%    none

% Lop off S cone direction if we are just doing L and M
if theDimension == 2 & size(directionCoding,1) > 2
   directionCoding(3:end,:) = []; 
end

% Apply the max contrast to each direction coding column
maxContDir   = bsxfun(@times,directionCoding,maxContrastPerDir);

% Make matricies the the size of the total desired stim.
fullContDir  = repelem(maxContDir,1,length(contrastSpacing));
fullContCode = repmat(contrastSpacing,1,length(maxContrastPerDir));

% apply the constrast spacining to the constast directions 
stimCombinations = bsxfun(@times,fullContDir,fullContCode);

end