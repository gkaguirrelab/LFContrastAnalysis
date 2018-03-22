%% analyzeLFContrast
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Convenience variables
projectName = 'LFContrastAnalysis';
flywheelName = 'LFContrast';

%% Analysis labels that we are going to go and get
fmriprepLabel = 'fmriprep 02/09/2018 11:40:55';
neuropythyLabel = 'retinotopy-templates 11/22/2017 13:21:46';
[fwInfo] = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,fmriprepDir, 'verbose', true, 'searchDir', projectDir, 'nodownload', true);


%% Set up params.
subjID       = 'sub-HEROgka1';
session      = 'ses-201709191435';
funcRuns     = {'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

%% Apply Warping to MNI space            
% Set up vars in order to run applyANTsWarpToData

% nifti input volumes (can be output of benson atlas or any nifti needed to
% be tranformed into preproc space) 
inRetFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};

% path to the retinotopy files

path2input   = ['~/Documents/flywheel/retAtlas/',subjID];
path2input = fullfile(getpref(projectName,'projectRootDir'),'neuropythy');

path2ref     = ['~/Documents/flywheel/fmriprep/',[subjID '_' ,'/',session,'/func'];
refFileName  = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = ['~/Documents/flywheel/fmriprep/', subjID, '/', session, '/anat'];
warpFileName = 'sub-HEROgka1_ses-201709191435_T1w_target-MNI152NLin2009cAsym_warp.h5';

% load ecc nifti file
eccenPos = find(~cellfun(@isempty,strfind(inRetFiles,'eccen')));
[~,tempName,~] = fileparts(inRetFiles{eccenPos});
[~,outName,~] = fileparts(tempName);
eccenFileName = fullfile(path2input,[outName '.nii.gz']);
eccen = MRIread(eccenFileName);

% load areas nifti file
areasPos = find(~cellfun(@isempty,strfind(inRetFiles,'areas')));
[~,tempName,~] = fileparts(inRetFiles{areasPos});
[~,outName,~] = fileparts(tempName);
areasFileName = fullfile(path2input,[outName,'.nii.gz']);
areas = MRIread(areasFileName);

% could add polar angle here but required current analysis
areaNum = 1;
eccenRange = [3 20];
[~,maskSaveName] = makeMaskFromRetino(eccen,areas,areaNum,eccenRange,path2input);

files2warp = {'HERO_gka1_T1.nii.gz',maskSaveName};

for ii = 1:length(files2warp)
    % input file
    inFile = fullfile(path2input,files2warp{ii});
    
    % output file
    [~,tempName,~] = fileparts(inFile);
    [~,outName,~] = fileparts(tempName);
    outFile = fullfile(path2input,[outName '_MNI_resampled.nii.gz']);
    
    % reference file
    refFile = fullfile(path2ref,refFileName);
    
    
    % warp file
    warpFile = fullfile(path2warp,warpFileName);
    if ~exist(outFile)
        applyANTsWarpToData(inFile, outFile, warpFile, refFile);
    end
end

%% Extract Signal from voxels
% Load mask nifti
[~,tmpName,~] = fileparts(maskSaveName);
[~,outName,~] = fileparts(tmpName);
maskOutFileName = fullfile(path2input,[outName '_MNI_resampled.nii.gz']);
mask = MRIread(maskOutFileName);
maskVol = mask.vol;

% make full file path to funvtional runs
funcRuns = fullfile(path2ref,funcRuns);

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