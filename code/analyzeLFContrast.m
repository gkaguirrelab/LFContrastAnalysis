%% Analyze LFContrast Data.
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Set up params.
subjID       = 'sub-HEROgka1';
session      = 'ses-201709191435';
funcRuns     = {'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};
% Run applyRetAtlas2Functional

% retinotopy file inputs
inRetFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz'};
% path to the retinotopy files
path2input   = ['~/Documents/flywheel/retAtlas/',subjID];
path2ref     = ['~/Documents/flywheel/fmriprep/',subjID,'/',session,'/func'];
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = ['~/Documents/flywheel/fmriprep/', subjID, '/', session, '/anat'];
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_space-MNI152NLin2009cAsym_warp.h5';

for ii = 1:length(inRetFiles)
    % input file
    inFile = fullfile(path2input,inRetFiles{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(path2input,[outName '_MNI_resampled.nii.gz']);
    
    % reference file
    refFile = fullfile(path2ref,refFileName);
    
    % warp file
    warpFile = fullfile(path2warp,warpFileName);
    if ~exist(outFile)
        applyANTsWarpToData(inFile, outFile, warpFile, refFile)
    end
end

%% Extract Signal from voxels
% load ecc
eccenPos = find(~cellfun(@isempty,strfind(inRetFiles,'eccen')));
[~,tempName,~] = fileparts(inRetFiles{eccenPos});
[~,outName,~] = fileparts(tempName);
eccenOutFileName = fullfile(path2input,[outName '_MNI_resampled.nii.gz'])
eccen = MRIread(eccenOutFileName);
% load areas
areasPos = find(~cellfun(@isempty,strfind(inRetFiles,'areas')));
[~,tempName,~] = fileparts(inRetFiles{areasPos});
[~,outName,~] = fileparts(tempName);
areasOutFileName = fullfile(path2input,[outName '_MNI_resampled.nii.gz'])
areas = MRIread(areasOutFileName);

% get maps
eccenMap = eccen.vol;
areasMap = areas.vol;

areaVal   = 1; % 1 = v1 2 = v2 3 = v3
eccenThresh = 30;
funcRuns = strcat(path2ref,funcRuns);
meanSignal = extractMeanSignalFromROI(funcRuns,areasMap,eccenMap,areaVal, eccenThresh)

