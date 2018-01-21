function path2ResampAtlas = applyRetAtlas2Functional(retFiles,warpFile,params,varargin)
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

p = inputParser;
p.addParameter('verbose',1,@isnumeric);
p.addParameter('dimensions',3,@isnumeric);

for ii = 1:length(retFiles)
p.parse(varargin{:});/home/mbarnett/Documents/flywheel/retAtlas/sub-HEROgka1/HERO_gka1_native.template_areas_MNIresampled2func.nii.gz 

cmd = antsApplyTransforms -d 3 -o /home/mbarnett/Documents/flywheel/retAtlas/sub-HEROgka1/HERO_gka1_native.template_areas_MNIresampled2func.nii.gz -v 1 -t /home/mbarnett/Documents/flywheel/fmriprep/sub-HEROgka1/ses-201709191435/anat/sub-HEROgka1_ses-201709191435_T1w_target-MNI152NLin2009cAsym_warp.h5 -i /home/mbarnett/Documents/flywheel/retAtlas/sub-HEROgka1/HERO_gka1_native.template_areas.nii.gz -r /home/mbarnett/Documents/flywheel/fmriprep/sub-HEROgka1/ses-201709191435/func/sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz




end


