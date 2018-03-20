% warpRetAtlasToFunc.m
%
% Example of how to call applyRetAtlas2Functional.m to transform the 
% retinotopy files from the benson atlas into MNI space using ANTs and the 
% warp file from fmriprep 
% 
% This can be deleted or turned into a demo -- not needed for analyzeLFContrast

% retinotopy file inputs -- from the benson atlas gear on flywheel 
inFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz'};
% path to the retinotopy files
path2input   = '~/Documents/flywheel/retAtlas/sub-HEROgka1';
path2ref     = '~/Documents/flywheel/fmriprep/sub-HEROgka1/ses-201709191435/func';
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = '~/Documents/flywheel/fmriprep/sub-HEROgka1/ses-201709191435/anat';
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_space-MNI152NLin2009cAsym_warp.h5';

for ii = 1:length(inFiles)
    % input file
    inFile = fullfile(path2input,inFiles{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(path2input,[outName '_MNI_resampled.nii.gz'])
    
    % reference file
    refFile = fullfile(path2ref,refFileName);
    
    % warp file
    warpFile = fullfile(path2warp,warpFileName);
    
    % run the warp
    applyANTsWarpToData(inFile, outFile, warpFile, refFile)
end