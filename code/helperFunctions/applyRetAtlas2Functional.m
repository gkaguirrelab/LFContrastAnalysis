function applyRetAtlas2Functional(retFiles,funcData,varargin)
% applyRetAtlas2Functional -- Apply the output of the Benson retinotopy
%                             gear to the functional data.
%
%   Details: This function uses the freesurfer command line tools to
%            register and resample the anatomical retinotopy template(s)
%            to match the functional data space. *Needs Freesurfer*
%
%   Inputs:
%       funcData            =  full path to the function run
%
%   varargin:
%
%   Output:
%       NONE
%
%   Call:
%       applyRetAtlas2Functional(retFiles,funcData,varargin)
%

% mab 2017 -- created

p = inputParser;

% Optional anaysis params
p.addParameter();
p.addParameter(,3,@isnumeric);
p.addParameter(,15,@isnumeric);
% Optional display params
p.addParameter('verbose',false,@islogical);
% parse
p.parse(varargin{:})


subjId = 'HERO_gka1'
% Set enviroment variable for the subject direcetory for freesurfer to
% locate revelant files

% this should changed to an input and smartly done with the parser once data directory
% structure is figured out in /tmp/flywheel/
pathToSubjDir = '/tmp/flywheel/HERO_gka1/flywheel/v0/output/fmriprep_output/freesurfer/';
setenv('SUBJECTS_DIR',pathToSubjDir)

% Create a registration file for the anatomical to the functional
funcData ='/tmp/flywheel/HERO_gka1/flywheel/v0/output/fmriprep_output/fmriprep/sub-HEROgka1/ses-201709191435/func/sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz';
subjId   = 'HEROgka1';
bbregName= 'register.dat';


system(['bbregister --mov ' funcData ' --bold --s ' subjId ' --init-fsl --reg ' fullFileReg]);

%% Set the file names
areaInName          =[subjId '_native.template_areas.nii.gz'];
areaOutName         = [subjId, '_reg2func.areas.nii.gz'];


%% Project area template to functional space
areaInFile                  = fullfile(params.sessionDir,'anat_templates',areaInName);

areas                 = fullfile(params.sessionDir,funcData{params.runNum},areaOutName);
cmd = ['mri_vol2vol --mov ' fullfile(params.sessionDir,funcData{params.runNum},[params.func '.nii.gz']) ...
    ' --targ ' areaInFile ' --o ' areas ...
    ' --reg ' bbregFile ' --inv --nearest'];
unix(cmd);
