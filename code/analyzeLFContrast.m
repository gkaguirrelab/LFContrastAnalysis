%% Analyze LFContrast Data.
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Set up params.
subjID       = 'HEROgka1';
session      = 'ses-201709191435';
funcRuns     = {'
% Run applyRetAtlas2Functional

% retinotopy file inputs
inFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz'};
% path to the retinotopy files
path2input   = ['~/Documents/flywheel/retAtlas/',subjID];
path2ref     = ['~/Documents/flywheel/fmriprep/',subjID,'/',session,'/func'];
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = ['~/Documents/flywheel/fmriprep/', subjID, '/', session, '/anat'];
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_space-MNI152NLin2009cAsym_warp.h5';

for ii = 1:length(inFiles)
    % input file
    inFile = fullfile(path2input,inFiles{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(path2input,[outName '_MNI_resampled.nii.gz'])
    
    % reference file
    refFile = fullfile(path2ref,refFileName);
    
    % warp file
    warpFile = fullfile(path2warp,warpFileName);
    
    applyANTsWarpToData(inFile, outFile, warpFile, refFile)
end

%% Extract Signal from voxels
% load ecc
eccenPos = find(~cellfun(@isempty,strfind(inFiles,'eccen')));
[~,tempName,~] = fileparts(inFiles{eccenPos});
[~,outName,~] = fileparts(tempName);
eccenOutFileName = fullfile(path2input,[outName '_MNI_resampled.nii.gz'])
eccen = MRIread(eccenOutFileName);
% load areas
areasPos = find(~cellfun(@isempty,strfind(inFiles,'areas')));
[~,tempName,~] = fileparts(inFiles{areasPos});
[~,outName,~] = fileparts(tempName);
areasOutFileName = fullfile(path2input,[outName '_MNI_resampled.nii.gz'])
areas = MRIread(areasOutFileName);

% get maps
eccenMap = eccen.vol;
areasMap = areas.vol;

areaVal   = 1; % 1 = v1 2 = v2 3 = v3
eccenThresh = 30;

for 
meanSignal = extractMeanSignalFromROI(timeSeries,areasMap,eccenMap,areaVal, eccenThresh)

%% Plot the data
figure; hold on
