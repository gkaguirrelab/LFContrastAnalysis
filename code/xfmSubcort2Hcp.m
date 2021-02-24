function [] = xfmSubcort2Hcp(origMgzPath, thalamSegPath, mni_2mm, empty_cifti_cortex_left, empty_cifti_cortex_right, template_cifti, TR, tmpDir, ciftiOutfile)
% This function takes the thalamic nuclei segmentation output of the 
%
% Inputs:
%   origMgzPath       = Unzipped folder containing the main recon-all results
%   thalamSegPath     = Label for creating subcortex
%   mni_2mm            = MNI template 2mm
%   empty_ciftify_cortex_left  = left cortex with all zeros
%   empty_ciftify_cortex_right = right cortex with all zeros
%   workdir           = workdir where intermediate files will be saved
%   outputdir         = outputdir where the final map will be saved
%   TR                = which TR to use (should be '1' [string])
%   tmpDir            = path to a temp folder to save intermediate files
%   template_cifti    = template cifti file 
%
% Outputs:
%   Saves 'subcortex_cifti.dtseries.nii' to outputdir 
%
% *NOTE*: MAB -- To run this in matlab you have to add the freefurfer
%         binaries to the path by running the following each time you start 
%         matlab:
%
%    bash_path=getenv ('PATH');
%    setenv( 'PATH',[bash_path,':/Applications/freesurfer7',':/Applications/freesurfer7/bin']);
%    setenv( 'PATH',[bash_path,':/Users/mbarnett/Documents/Code/ANTs/bin/ANTs-build/Examples/']);
%    setenv( 'PATH',[bash_path,':/Applications/workbench/bin_macosx64/']);
%    setenv('FREESURFER_HOME','/Applications/freesurfer7')
%  
%         Change '/Applications/freesurfer7' to match where you have 
%         freesurfer 
%
%         Change '/Users/mbarnett/Documents/Code/ANTs/bin/ANTs-build/Examples/'
%         to match where you have the ANTs binaries 
%
%         Change '/Applications/workbench/bin_macosx64/' to match where you 
%         have the workbench binaries 
%
% 
% OT  2/21    - wrote it in python
% MAB 2/24/21 - ported to matlab

%% Run MRI convert
% set up system command sting
% orig.mgz
orig_file_mgz = fullfile(origMgzPath, 'orig.mgz');
orig_file_nii = fullfile(tmpDir, 'orig.nii.gz');
mriConv_orig = ['mri_convert ' orig_file_mgz ' ' orig_file_nii];

% ThalamicNuclei mgz
segmentation_mgz   = fullfile(thalamSegPath, 'ThalamicNuclei.v12.T1.FSvoxelSpace.mgz');
segmentation_nifti = fullfile(tmpDir, 'ThalamicNuclei.v12.T1.FSvoxelSpace.nii.gz');
mriConv_thal = ['mri_convert ' segmentation_mgz ' ' segmentation_nifti];

% execute system calls
system(mriConv_orig)
system(mriConv_thal)

%% Run ANTS registration
% set up system command sting
antsReg = ['antsRegistration --verbose 1 --dimensionality 3 --float 0  ' ...
    '--collapse-output-transforms 1 --output [ ' tmpDir '/orig2MNI,' tmpDir '/orig2MNIWarped.nii.gz, ' tmpDir '/orig2MNIInverseWarped.nii.gz ] ' ...
    '--interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] ' ...
    '--initial-moving-transform [' mni_2mm ' , ' orig_file_nii ',1]' ...
    ' --transform Rigid[ 0.1 ] --metric MI[' mni_2mm ' , ' orig_file_nii ',1,32,Regular,0.25] '...
    '--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ' ...
    '--transform Affine[ 0.1 ] --metric MI[' mni_2mm  ' , ' orig_file_nii ',1,32,Regular,0.25] ' ...
    '--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox ' ...
    '--transform SyN[0.1,3,0] --metric CC[' mni_2mm ' , ' orig_file_nii ',1,4] ' ...
    '--convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox']

% execute system call
system(antsReg)

%% Apply warp
% set up system command sting
generic_affine = fullfile(tmpDir, 'orig2MNI0GenericAffine.mat')
warp = fullfile(tmpDir, 'orig2MNI1Warp.nii.gz')
segmentations_in_mni = fullfile(tmpDir, 'MNI_ThalamicNuclei.v12.T1.FSvoxelSpace.nii.gz')
antsApply = ['antsApplyTransforms -d 3 -i ' segmentation_nifti ' -r ' mni_2mm ' -t ' warp ' -t ' generic_affine ' -n NearestNeighbor -o ' segmentations_in_mni];

% execute system call
system(antsApply)

%% Create the CIFTI file
% set up system command sting
wb_cifti = ['wb_command -cifti-create-dense-from-template ' template_cifti ' ' ' ' ciftiOutfile '  -series ' TR ' 0 -volume-all ' segmentations_in_mni ...
           ' -metric CORTEX_LEFT ' empty_cifti_cortex_left ' -metric CORTEX_RIGHT ' empty_cifti_cortex_right];

% execute system call
system(wb_cifti)