%% UNCOMMENT AND RUN THIS IF YOU HAVE NOT DONE SO SINCE STARTIN MATLAB!!
%    bash_path=getenv ('PATH');
%    setenv( 'PATH',[bash_path,':/Applications/freesurfer7',':/Applications/freesurfer7/bin']);
%    setenv( 'PATH',[bash_path,':/Users/mbarnett/Documents/Code/ANTs/bin/ANTs-build/Examples/']);
%    setenv( 'PATH',[bash_path,':/Applications/workbench/bin_macosx64/']);
%    setenv('FREESURFER_HOME','/Applications/freesurfer7')

%% The good stuff
% subject ID
subjId = 'LZ23';

% Common path for all folders
baseDir = '/Users/mbarnett/labDropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/';

% Set up inputs to xfmSubcort2Hcp 
origMgzPath = fullfile(baseDir,'xfmToHCP',subjId); 
thalamSegPath =fullfile(baseDir,'xfmToHCP',subjId); 
mni_2mm  = fullfile(baseDir,'xfmToHCP','templates','MNI152_T1_2mm.nii.gz'); 
empty_cifti_cortex_left = fullfile(baseDir,'xfmToHCP','templates','empty_left.func.gii');
empty_cifti_cortex_right = fullfile(baseDir,'xfmToHCP','templates','empty_right.func.gii');
template_cifti = fullfile(baseDir,'xfmToHCP','templates','template_cifti.nii');
TR = '1';
tmpDir = fullfile(baseDir,'xfmToHCP','tmpDir'); 
ciftiOutfile = fullfile(baseDir,[subjId '_Thalamic_Seg_cifti.dtseries.nii']);

% Run the tranfromation from idividual subject to cifti space 
xfmSubcort2Hcp(origMgzPath, thalamSegPath, mni_2mm, empty_cifti_cortex_left, empty_cifti_cortex_right, template_cifti, TR, tmpDir, ciftiOutfile)