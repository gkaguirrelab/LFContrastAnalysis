% Run applyRetAtlas2Functional

% retinotopy file inputs
inFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz'};
% path to the retinotopy files
path2input   = '/home/mbarnett/Documents/flywheel/retAtlas/sub-HEROgka1';
path2ref     = '/home/mbarnett/Documents/flywheel/fmriprep/sub-HEROgka1/ses-201709191435/func';
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = '/home/mbarnett/Documents/flywheel/fmriprep/sub-HEROgka1/ses-201709191435/anat';
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_target-MNI152NLin2009cAsym_warp.h5';

for ii = 1:length(inFiles)
    % input file
    inFile = fullfile(path2input,inFiles{ii})
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(path2ret,[outName '_MNI_resampled.nii.gz']);
    
    % reference file
    refFile = fullfile(path2ref,refFileName);
    
    % warp file
    warpFile = fullfile(path2warp,warpFileName);
    
    [] = applyANTsWarpToData(inFile, outFile, warpFile, refFile)
end