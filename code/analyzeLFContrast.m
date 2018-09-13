% Convenience variables
analysisParams.projectName       = 'LFContrastAnalysis';
analysisParams.flywheelName      = 'LFContrast';
analysisParams.subjID            = 'sub-HEROGKA1';
analysisParams.expSubjID         = 'HERO_gka1';
analysisParams.session           = 'ses-ResearchAguirre';
analysisParams.sessionFolderName = 'HERO_GKA1_2018-07-28';
analysisParams.sessionDate       = '2018-07-28';
analysisParams.sessionNumber     = 'session_1';
analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
analysisParams.showPlots         = true;

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
analysisParams.areaNum     = 1;
analysisParams.eccenRange  = [3 20];

% Define the TR
analysisParams.TR = 0.800;
analysisParams.baselineCondNum = 6;
analysisParams.timeStep = 1/100;
analysisParams.generateIAMPPlots = true;


%Paramters for the QCM fit to IAMP:
analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
analysisParams.directionCoding = [1,1,1,0;-1,1,0,1;0,0,0,0]; %this 1 = L-M 2 = L+M 3 = L 4 = M;
analysisParams.maxContrastPerDir = [0.06,0.40,0.10,0.10]; % max contrast in the same order as above
analysisParams.theDimension = 2;

 %plotting params
 analysisParams.numSamples = 25;

[cleanRunData, analysisParams] = getTimeCourse(analysisParams);


[analysisParams,paramsQCMFit, meanIAMPBetas] = runIAMP_QCM(analysisParams,cleanRunData);

