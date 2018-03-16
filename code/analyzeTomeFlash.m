%%Analyze Flash Data.
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Set up params.
subjID       = 'sub-TOME3016';
session      = 'ses-Session2';
funcRuns     = {'sub-TOME3016_ses-Session2_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                'sub-TOME3016_ses-Session2_task-tfMRIFLASHPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

%% Apply Warping to MNI space            
% Set up vars in order to run applyANTsWarpToData

% nifti input volumes (can be output of benson atlas or any nifti needed to
% be tranformed into preproc space) 
inRetFiles = {'tome_3016_native.template_angle.nii.gz','tome_3016_native.template_areas.nii.gz','tome_3016_native.template_eccen.nii.gz',};

% path to the retinotopy files
path2input   = ['~/Documents/flywheel/retAtlas/',subjID];
path2ref     = ['~/Documents/flywheel/fmriprep/',subjID,'/',session,'/func'];
refFileName  = 'sub-TOME3016_ses-Session2_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
path2warp    = ['~/Documents/flywheel/fmriprep/', subjID, '/', session, '/anat'];
warpFileName = 'sub-TOME3016_ses-Session2_T1w_target-MNI152NLin2009cAsym_warp.h5';

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
eccenRange = [0 20];
[~,maskSaveName] = makeMaskFromRetino(eccen,areas,areaNum,eccenRange,path2input);

files2warp = {'tome_3016_T1.nii.gz',maskSaveName};

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

% get the mean signal in the mask for all time points and runs
meanSignal = extractMeanSignalFromMask(funcRuns,maskVol);

% convert to percent signal change relative to the mean of the run 
meanOfRun = repmat(mean(meanSignal),[size(meanSignal,1),1]);
PSC = 100*((meanSignal - meanOfRun)./meanOfRun);

figure; hold on
title('Time Course for TOME-3016') 
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
legend('Run 1','Run 2','Mean of Runs')
ylabel('Percent Signal Change (V1 0-20 deg)')
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