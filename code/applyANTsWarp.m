function [outFile] = applyANTsWarp(files2warp,refFileName,warpFileName)
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

for ii = 1:length(files2warp)
    % input file
    inFile = fullfile(retinoPath,files2warp{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);

    if ~exist(outFile)
        applyANTsWarpToData(inFile, outFile, warpFile, refFile);
    else
        [~,fileName,~] = fileparts(outFile);
        display(sprintf('%s already exist in the specified directory',fileName));
    end
end

end