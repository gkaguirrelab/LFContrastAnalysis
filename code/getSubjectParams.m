function analysisParams = getSubjectParams(subjId,varargin)
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('subjID',@ischar);
p.addParameter('roiType','V1',@ischar);
p.parse(subjId,varargin{:});

switch subjId
    case 'LZ23'
        % Convenience variables
        analysisParams.projectName       = 'LFContrastAnalysis';
        analysisParams.flywheelName      = 'LFContrast';
        analysisParams.subjID            = 'sub-LZ23';
        analysisParams.expSubjID         = 'LZ23';
        analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
        analysisParams.sessionFolderName = {'LZ23_2018-10-13','LZ23_2018-10-14'};
        analysisParams.sessionDate       = {'2018-10-13','2018-10-14'};
        analysisParams.sessionNumber     = {'session_1','session_1'};
        analysisParams.sessionNickname   = 'Original';
        analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
        analysisParams.showPlots         = true;
        analysisParams.projectNickname   = 'MRContrastResponseFunction';
        
        % Brain mask of function run for the reference volume in ANTs step
        analysisParams.refFileName  = 'sub-LZ23_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
        % output files of Neuropythy (retinotopy template)
        analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
        % warp file name (product of running fmriprep)
        analysisParams.warpFileName = 'sub-LZ23_T1w_target-MNI152NLin2009cAsym_warp.h5';
        
        % Paramters for the QCM fit to IAMP:
        analysisParams.theDimension = 2;
        analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
        analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5];
        analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM');
        analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
        analysisParams.numDirPerSession = 4;
        analysisParams.numRunsPerSession = 10;
        
        % Clip fisrt 2 TRs from time series?
        analysisParams.numClipFramesStart = 0;
        analysisParams.numClipFramesEnd   = 2;
        analysisParams.expLengthTR        = 360;
        
    case 'LZ23_replication'
        % Convenience variables
        analysisParams.projectName       = 'LFContrastAnalysis';
        analysisParams.flywheelName      = 'LFContrast';
        analysisParams.subjID            = 'sub-LZ23';
        analysisParams.expSubjID         = 'LZ23';
        analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
        analysisParams.sessionFolderName = {'LZ23_2019-03-10','LZ23_2019-03-09'};
        analysisParams.sessionDate       = {'2019-03-10','2019-03-09'};
        analysisParams.sessionNumber     = {'session_1','session_1'};
        analysisParams.sessionNickname   = 'Replication';
        analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
        analysisParams.showPlots         = true;
        analysisParams.projectNickname   = 'MRCRF';
        
        % Brain mask of function run for the reference volume in ANTs step
        analysisParams.refFileName  = 'sub-LZ23_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
        % output files of Neuropythy (retinotopy template)
        analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
        % warp file name (product of running fmriprep)
        analysisParams.warpFileName = 'sub-LZ23_T1w_target-MNI152NLin2009cAsym_warp.h5';
        
        % Paramters for the QCM fit to IAMP:
        analysisParams.theDimension = 2;
        analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
        analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5];
        analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM');
        analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
        analysisParams.numDirPerSession = 4;
        analysisParams.numRunsPerSession = 10;
        
        % Clip fisrt 2 TRs from time series?
        analysisParams.numClipFramesStart = 0;
        analysisParams.numClipFramesEnd   = 0;
        analysisParams.expLengthTR        = 360;
        
    case 'KAS25'
        % Convenience variables
        analysisParams.projectName       = 'LFContrastAnalysis';
        analysisParams.flywheelName      = 'LFContrast';
        analysisParams.subjID            = 'sub-KAS25';
        analysisParams.expSubjID         = 'KAS25';
        analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
        analysisParams.sessionFolderName = {'KAS25_2018-10-13','KAS25_2018-10-20'};
        analysisParams.sessionDate       = {'2018-10-13','2018-10-20'};
        analysisParams.sessionNumber     = {'session_1','session_1'};
        analysisParams.sessionNickname   = 'Original';
        analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
        analysisParams.showPlots         = true;
        analysisParams.projectNickname       = 'MRContrastResponseFunction';
        
        % Brain mask of function run fo r the reference volume in ANTs step
        analysisParams.refFileName  = 'sub-KAS25_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
        % output files of Neuropythy (retinotopy template)
        analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
        % warp file name (product of running fmriprep)
        analysisParams.warpFileName = 'sub-KAS25_T1w_target-MNI152NLin2009cAsym_warp.h5';
        
        % Paramters for the QCM fit to IAMP:
        analysisParams.theDimension = 2;
        analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
        analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5];
        analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM');
        analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
        analysisParams.numDirPerSession = 4;
        
        % Clip fisrt 2 TRs from time series?
        analysisParams.numClipFramesStart = 0;
        analysisParams.numClipFramesEnd   = 2;
        analysisParams.expLengthTR        = 360;
        
    case 'KAS25_replication'
        % Convenience variables
        analysisParams.projectName       = 'LFContrastAnalysis';
        analysisParams.flywheelName      = 'LFContrast';
        analysisParams.subjID            = 'sub-KAS25';
        analysisParams.expSubjID         = 'KAS25';
        analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
        analysisParams.sessionFolderName = {'KAS25_2019-03-16','KAS25_2019-03-17'};
        analysisParams.sessionDate       = {'2019-03-16','2019-03-17'};
        analysisParams.sessionNumber     = {'session_1','session_1'};
        analysisParams.sessionNickname   = 'Replication';
        analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
        analysisParams.showPlots         = true;
        analysisParams.projectNickname   = 'MRCRF';
        
        % Brain mask of function run fo r the reference volume in ANTs step
        analysisParams.refFileName  = 'sub-KAS25_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
        % output files of Neuropythy (retinotopy template)
        analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
        % warp file name (product of running fmriprep)
        analysisParams.warpFileName = 'sub-KAS25_T1w_target-MNI152NLin2009cAsym_warp.h5';
        
        % Paramters for the QCM fit to IAMP:
        analysisParams.theDimension = 2;
        analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
        analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5];
        analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM');
        analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
        analysisParams.numDirPerSession = 4;
        
        % Clip fisrt 2 TRs from time series?
        analysisParams.numClipFramesStart = 0;
        analysisParams.numClipFramesEnd   = 0;
        analysisParams.expLengthTR        = 360;
        
    case 'AP26'
        % Convenience variables
        analysisParams.projectName       = 'LFContrastAnalysis';
        analysisParams.flywheelName      = 'LFContrast';
        analysisParams.subjID            = 'sub-AP26';
        analysisParams.expSubjID         = 'AP26';
        analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
        analysisParams.sessionFolderName = {'AP26_2018-10-27','AP26_2018-10-21'};
        analysisParams.sessionDate       = {'2018-10-27','2018-10-21'};
        analysisParams.sessionNumber     = {'session_1','session_1'};
        analysisParams.sessionNickname   = 'Original';
        analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
        analysisParams.showPlots         = true;
        analysisParams.projectNickname   = 'MRContrastResponseFunction';
        
        % Brain mask of function run for the reference volume in ANTs step
        analysisParams.refFileName  = 'sub-AP26_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
        % output files of Neuropythy (retinotopy template)
        analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
        % warp file name (product of running fmriprep)
        analysisParams.warpFileName = 'sub-AP26_T1w_target-MNI152NLin2009cAsym_warp.h5';
        
        % Paramters for the QCM fit to IAMP:
        analysisParams.theDimension = 2;
        analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
        analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5];
        analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM');
        analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
        analysisParams.numDirPerSession = 4;
        
        % Clip fisrt 2 TRs from time series?
        analysisParams.numClipFramesStart = 0;
        analysisParams.numClipFramesEnd   = 2;
        analysisParams.expLengthTR        = 360;
        
    case 'AP26_replication'
        % Convenience variables
        analysisParams.projectName       = 'LFContrastAnalysis';
        analysisParams.flywheelName      = 'LFContrast';
        analysisParams.subjID            = 'sub-AP26';
        analysisParams.expSubjID         = 'AP26';
        analysisParams.session           = {'ses-ResearchAguirre','ses-ResearchAguirre'};
        analysisParams.sessionFolderName = {'AP26_2019-03-16','AP26_2019-03-31'};
        analysisParams.sessionDate       = {'2019-03-16','2019-03-31'};
        analysisParams.sessionNumber     = {'session_1','session_1'};
        analysisParams.sessionNickname   = 'Replication';
        analysisParams.sessionDir        = fullfile(getpref('LFContrastAnalysis','projectRootDir'),analysisParams.sessionFolderName);
        analysisParams.showPlots         = true;
        analysisParams.projectNickname   = 'MRCRF';
        
        % Brain mask of function run for the reference volume in ANTs step
        analysisParams.refFileName  = 'sub-AP26_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
        % output files of Neuropythy (retinotopy template)
        analysisParams.retinoFiles = {'rt_sub000_native.template_angle.nii.gz','rt_sub000_native.template_areas.nii.gz','rt_sub000_native.template_eccen.nii.gz',};
        % warp file name (product of running fmriprep)
        analysisParams.warpFileName = 'sub-AP26_T1w_target-MNI152NLin2009cAsym_warp.h5';
        
        % Paramters for the QCM fit to IAMP:
        analysisParams.theDimension = 2;
        analysisParams.contrastCoding = [1, .5, .25, .125, .0625];
        analysisParams.LMVectorAngles = [ -45, 45, 0, 90, -22.5, 22.5, 67.5, 112.5];
        analysisParams.directionCoding = vectorAngle2LMScontrast(analysisParams.LMVectorAngles,'LM');
        analysisParams.maxContrastPerDir = [0.12,0.60,0.14,0.22,0.085,0.20,0.40,0.13]; % max contrast in the same order as above
        analysisParams.numDirPerSession = 4;
        
        % Clip fisrt 2 TRs from time series?
        analysisParams.numClipFramesStart = 0;
        analysisParams.numClipFramesEnd   = 0;
        analysisParams.expLengthTR        = 360;
        
end


%% Add common params

switch p.Results.roiType
    case 'V1'
        % Info needed to make the V1 mask  from benson maps
        analysisParams.areaNum     = '1';
        analysisParams.eccenRange  = [0 20];
        analysisParams.anglesRange  = [0 180];
        analysisParams.hemisphere   = 'combined';
        analysisParams.threshold    = 0.9;
        analysisParams.useROI = false;
        
    case 'V2'
        analysisParams.areaNum     = '2';
        analysisParams.eccenRange  = [0 20];
        analysisParams.anglesRange  = [0 180];
        analysisParams.hemisphere   = 'combined';
        analysisParams.threshold    = 0.9;
        analysisParams.useROI = false;
        
    case 'V3'
        analysisParams.areaNum     = '3';
        analysisParams.eccenRange  = [0 20];
        analysisParams.anglesRange  = [0 180];
        analysisParams.hemisphere   = 'combined';
        analysisParams.threshold    = 0.9;
        analysisParams.useROI = false;
    case 'EVC'
        analysisParams.useROI = true;
        % rois located in Dropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/KAS25_LGN_ROI.dscalar.nii
        analysisParams.ROI    = ['combined.EVC.dscalar.nii'];    
    case 'LGN'
        analysisParams.useROI = true;
        % rois located in Dropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/KAS25_LGN_ROI.dscalar.nii
        analysisParams.ROI    = [analysisParams.expSubjID '_LGN_ROI.dscalar.nii'];
    case 'V4'
        analysisParams.useROI = true;
        % rois located in Dropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/KAS25_LGN_ROI.dscalar.nii
        analysisParams.ROI    = [analysisParams.expSubjID '_wang_V4.dscalar.nii'];
    case 'VO1'
        analysisParams.useROI = true;
        % rois located in Dropbox/MELA_analysis/LFContrastAnalysis/MNI_ROIs/KAS25_LGN_ROI.dscalar.nii
        analysisParams.ROI    = [analysisParams.expSubjID '_wang_VO1.dscalar.nii'];
end


% Define the TR
analysisParams.TR = 0.800;
analysisParams.baselineCondNum = 6;
analysisParams.timeStep = 1/100;
analysisParams.generateIAMPPlots = false;
analysisParams.generateCrossValPlots = false;
analysisParams.blockDuration = 12; %seconds
analysisParams.numFramesPerBlock = analysisParams.TR * analysisParams.blockDuration;

% Plotting params
analysisParams.numSamples = 25;

% create timebase
totalTime =analysisParams.expLengthTR*analysisParams.TR;
deltaT = analysisParams.TR;
analysisParams.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);

