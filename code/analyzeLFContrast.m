% % LZ23 
% Convenience variables
analysisParams.projectName       = 'LFContrastAnalysis';
analysisParams.flywheelName      = 'LFContrast';
analysisParams.subjID            = 'sub-LZ23';
analysisParams.expSubjID         = 'LZ23';
analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
analysisParams.sessionFolderName = {'LZ23_2018-10-13','LZ23_2018-10-14'};
analysisParams.sessionDate       = {'2018-10-13','2018-10-14'};
analysisParams.sessionNumber     = {'session_1','session_1'};
analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
analysisParams.showPlots         = true;

% Brain mask of function run for the reference volume in ANTs step
analysisParams.refFileName  = 'sub-LZ23_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
% output files of Neuropythy (retinotopy template)
analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
% warp file name (product of running fmriprep)
analysisParams.warpFileName = 'sub-LZ23_T1w_target-MNI152NLin2009cAsym_warp.h5';

% Paramters for the QCM fit to IAMP:
analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5]; 
analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM'); 
analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
analysisParams.theDimension = 2;

% Clip fisrt 2 TRs from time series?
% if no clipping then put 0;
analysisParams.numClipFramesStart = 0;
analysisParams.numClipFramesEnd   = 2;

% Make mask from the area and eccentricity maps
analysisParams.areaNum     = 1;
analysisParams.eccenRange  = [21 26];

% Define the TR
analysisParams.TR = 0.800;
analysisParams.baselineCondNum = 6;
analysisParams.timeStep = 1/100;
analysisParams.generateIAMPPlots = false;

% Plotting params
 analysisParams.numSamples = 25;

% Get the cleaned time series
[fullCleanData, analysisParams] = getTimeCourse(analysisParams);

% Run the IAMP/QCM model
[analysisParams,paramsQCMFit, meanIAMPBetas, semIAMPBetas,packetPocket,paramsFitIAMP,fitResponseStructQCM, baselineBetas] = runIAMP_QCM(analysisParams,fullCleanData);

% Plot the CRF from the IAMP and QCM fits
plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit);

% Plot isoresponce contour
thresholds = [.1, .2, .4];
hdl = plotIsoresponse(analysisParams,meanIAMPBetas,paramsQCMFit,thresholds);

% Use QCM fit to IAMP to predict timecourse.
plotQCMtimecourse(paramsFitIAMP,packetPocket,meanIAMPBetas,analysisParams,fitResponseStructQCM,baselineBetas);


%% STORED VARS FOR OTHER SUBJECTS:

%% GKA1
%Convenience variables
% analysisParams.projectName       = 'LFContrastAnalysis';
% analysisParams.flywheelName      = 'LFContrast';
% analysisParams.subjID            = 'sub-HEROGKA1';
% analysisParams.expSubjID         = 'HERO_gka1';
% analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
% analysisParams.sessionFolderName = {'HERO_GKA1_2018-07-28','HERO_GKA1_2018-08-21'};
% analysisParams.sessionDate       = {'2018-07-28','2018-08-21'};
% analysisParams.sessionNumber     = {'session_1','session_1'};
% analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
% analysisParams.showPlots         = true;
% analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
% analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5]; 
% analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM'); 
% analysisParams.maxContrastPerDir = [0.06,0.40,0.10,0.10,0.085,0.24,0.20,0.10]; % max contrast in the same order as above



% % LZ23 
% % Convenience variables
% analysisParams.projectName       = 'LFContrastAnalysis';
% analysisParams.flywheelName      = 'LFContrast';
% analysisParams.subjID            = 'sub-LZ23';
% analysisParams.expSubjID         = 'LZ23';
% analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
% analysisParams.sessionFolderName = {'LZ23_2018-10-13','LZ23_2018-10-14'};
% analysisParams.sessionDate       = {'2018-10-13','2018-10-14'};
% analysisParams.sessionNumber     = {'session_1','session_1'};
% analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
% analysisParams.showPlots         = true;
%
% % Brain mask of function run for the reference volume in ANTs step
% analysisParams.refFileName  = 'sub-LZ23_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
% % output files of Neuropythy (retinotopy template)
% analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
% % warp file name (product of running fmriprep)
% analysisParams.warpFileName = 'sub-LZ23_T1w_target-MNI152NLin2009cAsym_warp.h5';
%
% % Paramters for the QCM fit to IAMP:
% analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
% analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5]; 
% analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM'); 
% analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
% analysisParams.theDimension = 2;


% % KAS25
% % Convenience variables
% analysisParams.projectName       = 'LFContrastAnalysis';
% analysisParams.flywheelName      = 'LFContrast';
% analysisParams.subjID            = 'sub-KAS25';
% analysisParams.expSubjID         = 'KAS25';
% analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
% analysisParams.sessionFolderName = {'KAS25_2018-10-13','KAS25_2018-10-20'};
% analysisParams.sessionDate       = {'2018-10-13','2018-10-20'};
% analysisParams.sessionNumber     = {'session_1','session_1'};
% analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
% analysisParams.showPlots         = true;
% 
% % Brain mask of function run for the reference volume in ANTs step
% analysisParams.refFileName  = 'sub-KAS25_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
% % output files of Neuropythy (retinotopy template)
% analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
% % warp file name (product of running fmriprep)
% analysisParams.warpFileName = 'sub-KAS25_T1w_target-MNI152NLin2009cAsym_warp.h5';
% 
% % Paramters for the QCM fit to IAMP:
% analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
% analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5]; 
% analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM'); 
% analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
% analysisParams.theDimension = 2;




