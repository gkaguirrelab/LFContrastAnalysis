% Convenience variables
analysisParams.projectName       = 'LFContrastAnalysis';
analysisParams.flywheelName      = 'LFContrast';
analysisParams.subjID            = 'sub-HEROGKA1';
analysisParams.expSubjID         = 'HERO_gka1';
analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
analysisParams.sessionFolderName = {'HERO_GKA1_2018-07-28','HERO_GKA1_2018-08-21'};
analysisParams.sessionDate       = {'2018-07-28','2018-08-21'};
analysisParams.sessionNumber     = {'session_1','session_1'};
analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
analysisParams.showPlots         = true;

% Brain mask of function run for the reference volume in ANTs step
analysisParams.refFileName  = 'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
% output files of Neuropythy (retinotopy template)
analysisParams.retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};
% warp file name (product of running fmriprep)
analysisParams.warpFileName = 'sub-HEROGKA1_T1w_target-MNI152NLin2009cAsym_warp.h5';

% Clip fisrt 2 TRs from time series?
% if no clipping then put 0;
analysisParams.numClipFrames = 0;

% Make mask from the area and eccentricity maps
analysisParams.areaNum     = 1;
analysisParams.eccenRange  = [3 20];

% Define the TR
analysisParams.TR = 0.800;
analysisParams.baselineCondNum = 6;
analysisParams.timeStep = 1/100;
analysisParams.generateIAMPPlots = false;

% Paramters for the QCM fit to IAMP:
analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5]; 
analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM'); 
analysisParams.maxContrastPerDir = [0.06,0.40,0.10,0.10,0.085,0.24,0.20,0.10]; % max contrast in the same order as above
analysisParams.theDimension = 2;

% Plotting params
 analysisParams.numSamples = 25;

% Get the cleaned time series
[fullCleanData, analysisParams] = getTimeCourse(analysisParams);

% Run the IAMP/QCM model
[analysisParams,paramsQCMFit, meanIAMPBetas, semIAMPBetas,packetPocket,paramsFitIAMP,fitResponseStructQCM] = runIAMP_QCM(analysisParams,fullCleanData);

% Plot the CRF from the IAMP and QCM fits
plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit);

% Plot isoresponce contour
thresholds = [.25, .5, .75];
hdl = plotIsoresponse(analysisParams,meanIAMPBetas,paramsQCMFit,thresholds);

% Use QCM fit to IAMP to predict timecourse.
plotQCMtimecourse(paramsFitIAMP,packetPocket,meanIAMPBetas,analysisParams,fitResponseStructQCM);