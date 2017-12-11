function path2ResampAtlas = applyRetAtlas2Functional(retFiles,path2FuncData,params)
% applyRetAtlas2Functional -- Apply the output of the Benson retinotopy
%                             gear to the functional data.
%
%   Details: This function uses the freesurfer command line tools to
%            register and resample the anatomical retinotopy template(s)
%            to match the functional data space. *Needs Freesurfer*
%
%   Inputs:
%       funcData        =  full path to the function run
%       params          =  params.subjID -- subject ID as in the
%                                           freeserfuer directory 
%   varargin:
%
%   Output:
%       NONE
%
%   Call:
%       applyRetAtlas2Functional(retFiles,funcData,varargin)
%

% mab 2017 -- created

% Set enviroment variable for the subject direcetory for freesurfer to
% locate revelant files

% this should changed to an input and smartly done with the parser once data directory
% structure is figured out in /tmp/flywheel/

pathToSubjFS = ['/tmp/flywheel/' params.subjID '/flywheel/v0/output/fmriprep_output/freesurfer/'];
setenv('SUBJECTS_DIR',pathToSubjFS)

% Create a registration file for the anatomical to the functional
path2FuncData  = fullfile('tmp','flywheel', params.subjID, 'flywheel', 'v0', 'output','fmriprep_output','fmriprep',['sub-' params.subjID ], params.session, 'func');
funcFileName   = ['sub-' params.subjID '_' params.session '_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'];
funcFullFile   = [path2FuncData funcFileName];
fullFileReg    = ['/path/' params.regFileName]; % <-- NEED TO FIND A PLACE TO PUT THIS


% Option to run bbregister 
system(['bbregister --mov ' funcFullFile ' --bold --s ' subjId ' --init-fsl --reg ' params.regFileName  fullFileReg]);

%% Set the file names
areaInName          =[subjId '_native.template_areas.nii.gz'];
areaOutName         = [subjId, '_reg2func.areas.nii.gz'];


%% Project area template to functional space
areaInFile                  = fullfile(params.sessionDir,'anat_templates',areaInName);

areas                 = fullfile(params.sessionDir,path2FuncData{params.runNum},areaOutName);
cmd = ['mri_vol2vol --mov ' fullfile(params.sessionDir,path2FuncData{params.runNum},[params.func '.nii.gz']) ...
    ' --targ ' areaInFile ' --o ' areas ...
    ' --reg ' bbregFile ' --inv --nearest'];
unix(cmd);
