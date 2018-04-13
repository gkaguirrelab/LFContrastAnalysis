%%Analyze Flash Data.
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Set up params.
subjID         = 'sub-HEROgka1';
projectName    = 'LFContrastAnalysis';
flywheelName   = 'LFContrast';
session        = 'ses-0411181853PM';
functionalRuns = {'sub-HEROgka1_ses-0411181853PM_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

%% Analysis labels that we are going to go and get
fmriprepLabel   = 'fmriprep 04/12/2018 15:16:06';
neuropythyLabel = 'retinotopy-templates 11/22/2017 13:21:46';
fwInfo          = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,'', 'nodownload', true);
sessionDir      = fullfile(getpref('LFContrastAnalysis','projectRootDir'),[fwInfo.subject,'_', fwInfo.timestamp(1:10)]);

%% Apply Warping to MNI space            
% Set up vars in order to run applyANTsWarpToData

% nifti input volumes (can be output of benson atlas or any nifti needed to
% be tranformed into preproc space) 
inRetFiles = {'tome_3016_native.template_angle.nii.gz','tome_3016_native.template_areas.nii.gz','tome_3016_native.template_eccen.nii.gz',};

% path to the retinotopy files
path2input   = ['~/Documents/flywheel/retAtlas/',subjID];
path2ref     = ['~/Documents/flywheel /fmriprep/',subjID,'/',session,'/func'];
refFileName  = 'sub-TOME3016_ses-Session2_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = ['~/Documents/flywheel/fmriprep/', subjID, '/', session, '/anat'];
warpFileName = 'sub-TOME3016_ses-Session2_T1w_target-MNI152NLin2009cAsym_warp.h5';


numAcquisitions = length(functionalRuns);
% brain mask of function run for the reference volume in ANTs step
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';

% output files of Neuropythy (retinotopy template)
retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};

% warp file name (product of running fmriprep)
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_target-MNI152NLin2009cAsym_warp.h5';

% Set up paths to nifti and .h5 files
retinoPath     = fullfile(sessionDir,'neuropythy');
functionalPath = fullfile(sessionDir, 'fmriprep', subjID, session, 'func');
warpFilePath   = fullfile(sessionDir, 'fmriprep', subjID, session, 'anat');


%% Create restricted V1 mask

% load ecc nifti file
eccenPos       = find(~cellfun(@isempty,strfind(retinoFiles,'eccen')));
[~,tempName,~] = fileparts(retinoFiles{eccenPos});
[~,outName,~]  = fileparts(tempName);
eccenFileName  = fullfile(retinoPath,[outName '.nii.gz']);
eccen          = MRIread(eccenFileName);

% load areas nifti file
areasPos       = find(~cellfun(@isempty,strfind(retinoFiles,'areas')));
[~,tempName,~] = fileparts(retinoFiles{areasPos});
[~,outName,~]  = fileparts(tempName);
areasFileName  = fullfile(retinoPath,[outName,'.nii.gz']);
areas          = MRIread(areasFileName);

% make mask from the area and eccentricity maps
areaNum     = 1;
eccenRange  = [3 20];
[~,maskSaveName] = makeMaskFromRetino(eccen,areas,areaNum,eccenRange,retinoPath);

files2warp = {'HERO_gka1_T1.nii.gz',maskSaveName};
for ii = 1:length(files2warp)
    % input file
    inFile = fullfile(retinoPath,files2warp{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);
    
    % reference file
    refFile = fullfile(functionalPath,refFileName);
    
    % warp file
    warpFile = fullfile(warpFilePath,warpFileName);
    if ~exist(outFile)
        applyANTsWarpToData(inFile, outFile, warpFile, refFile);
    else
        [~,fileName,~] = fileparts(outFile);
        display(sprintf('%s already exist in the specified directory',fileName));
    end
end


%% Extract Signal from voxels
% Load mask nifti
maskPos       = find(~cellfun(@isempty,strfind(files2warp,'mask')));
[~,tempName,~] = fileparts(files2warp{maskPos});
[~,tmpName,~] = fileparts(maskSaveName);
[~,outName,~] = fileparts(tmpName);
maskOutFileName = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);
mask = MRIread(maskOutFileName);
maskVol = mask.vol;


% make full file path to funvtional runs
funcRuns = fullfile(functionalPath,functionalRuns);

% get the mean signal in the mask for all time points and runs
meanSignal = extractMeanSignalFromMask(funcRuns,maskVol);

% convert to percent signal change relative to the mean of the run 
meanOfRun = repmat(mean(meanSignal),[size(meanSignal,1),1]);
PSC = 100*((meanSignal - meanOfRun)./meanOfRun);

figure; hold on
title('Time Course for HEROgka1') 
timepoints = 0.8.*[1:size(PSC,1)]-0.8;
plot(timepoints,PSC,'--');
plot(timepoints,mean(PSC,2),'k', 'LineWidth',2);
minVal = -5;
maxVal = 5;
times = 0:12:(size(PSC,1).*0.8);
for ii = 1:length(times)
    plot([times(ii) times(ii)], [minVal maxVal],'Color',[0.5 0.5 0.5],'LineWidth',1);
end
ylim([-3 3]);
legend('Flash Run Data')
ylabel('Percent Signal Change (V1 3-20 deg)')
xlabel('Time (Seconds)')


%%plot stuff

% load the last param file (doesnt matter because all were the same
% load(dataParamFile);
% plotTimeCourse(meanSignal,block,responseStruct);

% %plot CRF
% figure;
% A = repmat(avgPerCond(end,:),[6,1]);
% B = 100*((avgPerCond - A)./A);
% plot([.8,.4,.2,.1,.05,0],mean(B,2))
% ylabel('Percent Signal Change')
% xlabel('Contrast Level')
% 
% title('Contrast Response Function')