%% analyzeLFContrast
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Convenience variables
projectName = 'LFContrastAnalysis';
flywheelName = 'LFContrast';
subjID       = 'sub-HEROgka1';
session      = 'ses-201709191435';

%% Analysis labels that we are going to go and get
fmriprepLabel = 'fmriprep 02/09/2018 11:40:55';
neuropythyLabel = 'retinotopy-templates 11/22/2017 13:21:46';
[fwInfo] = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,'', 'nodownload', true);
sessionDir = fullfile(getpref('LFContrastAnalysis','projectRootDir'),[fwInfo.subject,'_', fwInfo.timestamp(1:10)]);

%% Relevant Nifti names for analysis

% functional runs
funcRuns     = {'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};
            
% brain mask of function run for the reference volume in ANTs step
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';

% output files of Neuropythy (retinotopy template)
inRetFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};

% warp file name (product of running fmriprep) 
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_target-MNI152NLin2009cAsym_warp.h5';

%% Apply Warping to MNI space

% Set up vars in order to run applyANTsWarpToData
retinoPath     = fullfile(sessionDir,'neuropythy');
functionalPath = fullfile(sessionDir, 'fmriprep', subjID, session, 'func');
warpFilePath   = fullfile(sessionDir, 'fmriprep', subjID, session, 'anat');

% load ecc nifti file
eccenPos = find(~cellfun(@isempty,strfind(inRetFiles,'eccen')));
[~,tempName,~] = fileparts(inRetFiles{eccenPos});
[~,outName,~] = fileparts(tempName);
eccenFileName = fullfile(retinoPath,[outName '.nii.gz']);
eccen = MRIread(eccenFileName);

% load areas nifti file
areasPos = find(~cellfun(@isempty,strfind(inRetFiles,'areas')));
[~,tempName,~] = fileparts(inRetFiles{areasPos});
[~,outName,~] = fileparts(tempName);
areasFileName = fullfile(retinoPath,[outName,'.nii.gz']);
areas = MRIread(areasFileName);

% could add polar angle here but required current analysis
areaNum = 1;
eccenRange = [3 20];
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
    end
end

%% Extract Signal from voxels
% Load mask nifti
[~,tmpName,~] = fileparts(maskSaveName);
[~,outName,~] = fileparts(tmpName);
maskOutFileName = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);
mask = MRIread(maskOutFileName);
maskVol = mask.vol;

% make full file path to funvtional runs
funcRuns = fullfile(functionalPath,funcRuns);

% extract the mean signal from voxels
meanSignal = extractMeanSignalFromMask(funcRuns,maskVol);

% convert to percent signal change relative to the mean 
A = repmat(mean(meanSignal,1),[size(meanSignal,1),1]);
PSC = 100*((meanSignal - A)./A);


%% Get trial order info:
trialOrderDir = '~/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_TrialSequenceMR/MRContrastResponseFunction/DataFiles/HERO_gka1/2017-09-19/session_1';
trialOrderFiles = {'session_1_CRF_scan1.mat', 'session_1_scan2.mat', 'session_1_scan3.mat', 'session_1_scan4.mat', 'session_1_scan5.mat', 'session_1_scan6.mat'};

for jj = 1:length(trialOrderFiles)

    dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
    TR = 0.800;
    expParams = getExpParams(dataParamFile,TR);
    [avgPerCond(:,jj), blockOrder(:,jj)] = sortDataByConditions(PSC(:,jj),expParams);
    
end


%%plot stuff

% load the last param file (doesnt matter because all were the same
load(dataParamFile);
plotTimeCourse(meanSignal,block,responseStruct);

%plot CRF
figure; hold on
zeroCondMean = repmat(avgPerCond(6,:),[size(avgPerCond,2),1]);
CRF = avgPerCond -zeroCondMean; 
plot([.8,.4,.2,.1,.05,0],mean(CRF,2),'k','LineWidth',2)
plot([.8,.4,.2,.1,.05,0],CRF,'--o')
ylabel('Percent Signal Change (diff from zero cond)')
xlabel('Contrast Level')
legend('Mean of Runs', 'Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5', 'Run 6')
axis square
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Contrast Response Function')