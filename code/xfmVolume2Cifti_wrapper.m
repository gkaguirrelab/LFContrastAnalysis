%% UNCOMMENT AND RUN THIS IF YOU HAVE NOT DONE SO SINCE STARTIN MATLAB!!
%    bash_path=getenv ('PATH');
%    setenv( 'PATH',[bash_path,':/Applications/freesurfer7',':/Applications/freesurfer7/bin']);
%    setenv( 'PATH',[bash_path,':/Users/mbarnett/Documents/Code/ANTs/bin/ANTs-build/Examples']);
%    setenv( 'PATH',[bash_path,':/Applications/workbench/bin_macosx64']);
%    setenv('FREESURFER_HOME','/Applications/freesurfer7')


%% The good stuff

% Common path for all folders
baseDir = '/Users/mbarnett/labDropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/';

%% Set up inputs to xfmSubcort2Hcp 
% The movable nifti (infile)
refFileMovable = fullfile(baseDir,'ProbAtlas_v4','MNI152_T1_1mm.nii.gz'); 
niftiMovable =fullfile(baseDir,'ProbAtlas_v4','subj_vol_all','maxprob_vol_rh.nii.gz'); 

% The target nifti
nifitTagetVol  = fullfile(baseDir,'xfmToHCP','templates','MNI152_T1_2mm.nii.gz'); 

% templates and info of the cifti stuff
empty_cifti_cortex_left = fullfile(baseDir,'xfmToHCP','templates','empty_left.func.gii');
empty_cifti_cortex_right = fullfile(baseDir,'xfmToHCP','templates','empty_right.func.gii');
template_cifti = fullfile(baseDir,'xfmToHCP','templates','template_cifti.nii');
TR = '0';

% set a temp dir
tmpDir = fullfile(baseDir,'xfmToHCP','tmpDir'); 

% outfile name
ciftiOutfile = fullfile(baseDir,'maxprob_vol_rh.dtseries.nii');

% Run the tranfromation from idividual subject to cifti space 
xfmVolume2Cifti(refFileMovable, niftiMovable, nifitTagetVol, empty_cifti_cortex_left, empty_cifti_cortex_right, template_cifti, TR, tmpDir, ciftiOutfile)
