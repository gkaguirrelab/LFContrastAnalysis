function [] = applyANTsWarp()
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


%% Create restricted V1 mask
% load ecc nifti file
eccenPos       = find(~cellfun(@isempty,strfind(retinoFiles,'eccen')));
[~,tempName,~] = fileparts(retinoFiles{eccenPos});
[~,outName,~]  = fileparts(tempName);
eccenFileName  = fullfile(retinoPath,[outName '.nii.gz']);
eccen          = MRIread(eccenFileName);

% load areas nifti file
areasPos       = find(~cellfun(@isempty,strfind(retinoFiles,'areas')));
[~,tempName,~] = fileparts(retinoFiles{areasPos});
[~,outName,~]  = fileparts(tempName);
areasFileName  = fullfile(retinoPath,[outName,'.nii.gz']);
areas          = MRIread(areasFileName);

% make mask from the area and eccentricity maps
areaNum     = 1;
eccenRange  = [3 20];
[~,maskSaveName] = makeMaskFromRetino(eccen,areas,areaNum,eccenRange,retinoPath);

%% Apply the warp to the mask and T1 files using ANTs

files2warp = {'HERO_gka1_T1.nii.gz',maskSaveName};
for ii = 1:length(files2warp)
    % input file
    inFile = fullfile(retinoPath,files2warp{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);
    
    % reference file
    refFile = fullfile(functionalPath,refFileName);
    
    % warp file
    warpFile = fullfile(warpFilePath,warpFileName);
    if ~exist(outFile)
        applyANTsWarpToData(inFile, outFile, warpFile, refFile);
    else
        [~,fileName,~] = fileparts(outFile);
        display(sprintf('%s already exist in the specified directory',fileName));
    end
end