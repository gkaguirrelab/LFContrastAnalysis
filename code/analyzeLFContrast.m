 % Convenience variables
analysisParams.projectName       = 'LFContrastAnalysis';
analysisParams.flywheelName      = 'LFContrast';
analysisParams.subjID            = 'sub-HEROGKA1';
analysisParams.expSubjID         = 'HERO-gka1';
analysisParams.session           = 'ses-ResearchAguirre';
analysisParams.sessionFolderName = 'HERO_GKA1_2018-07-28';
analysisParams.sessionDate       = '2018-07-28';
analysisParams.sessionNumber     = 'session_1';
analysisParams.sessionDir = fullfile(getpref('LFContrastAnalysis','projectRootDir'),sessionFolderName);

% brain mask of function run for the reference volume in ANTs step
analysisParams.refFileName  = 'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
% output files of Neuropythy (retinotopy template)
analysisParams.retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};
% warp file name (product of running fmriprep)
analysisParams.warpFileName = 'sub-HEROGKA1_T1w_target-MNI152NLin2009cAsym_warp.h5';


% Clip fisrt 2 TRs from time series?
% if no clipping then put 0;
analysisParams.numClipFrames = 0;

% make mask from the area and eccentricity maps
analysisParams.areaNum     = 2;
analysisParams.eccenRange  = [3 20];

% Define the TR
analysisParams.TR = 0.800;

analysisParams = getTimeCourse(analysisParams);


%% Get trial order info:
trialOrderDir = fullfile(getpref(analysisParams.projectName,'melaDataPath'),HERO_gka1/2018-07-28/session_1');
trialOrderFiles = {'CRF_analysisParams.session_1_scan1.mat' ...
    'CRF_session_1_scan2.mat', ...
    'CRF_session_1_scan3.mat', ...
    'CRF_session_1_scan4.mat', ...
    'CRF_session_1_scan5.mat', ...
    'CRF_session_1_scan6.mat', ...
    'CRF_session_1_scan7.mat', ...
    'CRF_session_1_scan8.mat', ...
    'CRF_session_1_scan9.mat', ...
    'CRF_session_1_scan10.mat'};