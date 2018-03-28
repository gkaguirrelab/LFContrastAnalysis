%% analyzeLFContrast
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Convenience variables
projectName  = 'LFContrastAnalysis';
flywheelName = 'LFContrast';
subjID       = 'sub-HEROgka1';
session      = 'ses-201709191435';

%% Analysis labels that we are going to go and get
fmriprepLabel   = 'fmriprep 02/09/2018 11:40:55';
neuropythyLabel = 'retinotopy-templates 11/22/2017 13:21:46';
fwInfo          = getAnalysisFromFlywheel(flywheelName,fmriprepLabel,'', 'nodownload', true);
sessionDir      = fullfile(getpref('LFContrastAnalysis','projectRootDir'),[fwInfo.subject,'_', fwInfo.timestamp(1:10)]);

%% Relevant Nifti names for analysis

% functional runs
functionalRuns = {'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastAP_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

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

%% Apply the warp to the mask and T1 files using ANTs

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

% make full file path to functional runs
functionalRuns = fullfile(functionalPath,functionalRuns);

% extract the mean signal from voxels
[voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRuns,maskVol);
meanSignal = squeeze(mean(voxelTimeSeries,1));

% convert to percent signal change relative to the mean
meanMat = repmat(mean(meanSignal,1),[size(meanSignal,1),1]);
PSC = 100*((meanSignal - meanMat)./meanMat);


%% Get trial order info:
trialOrderDir = fullfile(getpref(projectName,'melaDataPath'),'/Experiments/OLApproach_TrialSequenceMR/MRContrastResponseFunction/DataFiles/HERO_gka1/2017-09-19/session_1');
trialOrderFiles = {'session_1_CRF_scan1.mat', 'session_1_scan2.mat', 'session_1_scan3.mat', 'session_1_scan4.mat', 'session_1_scan5.mat', 'session_1_scan6.mat'};

%% Construct the model object
temporalFit = tfeIAMP('verbosity','none');


%% Create a cell of stimulusStruct (one struct per run)
for jj = 1:length(trialOrderFiles)
    dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
    load(dataParamFile);
    TR = 0.800;
    expParams = getExpParams(dataParamFile,TR,'hrfOffset', false, 'stripInitialTRs', false);
    
    % make stimulus timebase
    totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
    deltaT = 800;
    stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
    
    % make stimulus values
    % Stim coding: 80% = 1, 40% = 2, 20% = 3, 10% = 4, 5% = 5, 0% = 6;
    stimulusStruct.values = zeros(length(unique(expParams(:,3))),totalTime/deltaT);
    for kk = 1:size(expParams,1)
        stimulusStruct.values(expParams(kk,3),expParams(kk,1):expParams(kk,2)) = 1;
    end
    % Set the number of instances. Right now hard-coded; should get this
    % from the stimulus information
    defaultParamsInfo.nInstances = 6;

    % Define HRF Copied from the t_BTRMBasic demo (a double gamma HRF)
    hrfParams.gamma1 = 6;   % positive gamma parameter (roughly, time-to-peak in secs)
    hrfParams.gamma2 = 12;  % negative gamma parameter (roughly, time-to-peak in secs)
    hrfParams.gammaScale = 10; % scaling factor between the positive and negative gamma componenets
    kernelStruct.timebase=stimulusStruct.timebase;
    % The timebase is converted to seconds within the function, as the gamma
    % parameters are defined in seconds.
    hrf = gampdf(kernelStruct.timebase/1000, hrfParams.gamma1, 1) - ...
        gampdf(kernelStruct.timebase/1000, hrfParams.gamma2, 1)/hrfParams.gammaScale;
    kernelStruct.values=hrf;
    % prepare this kernelStruct for use in convolution as a BOLD HRF
    kernelStruct.values=kernelStruct.values-kernelStruct.values(1);
    kernelStruct=normalizeKernelArea(kernelStruct);
        
    
    % make the stimulus portion of packet for fitting
    thePacket.stimulus = stimulusStruct;
    
    % add the response field
    responseStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
    thePacket.response.timebase = responseStruct.timebase;
    thePacket.response.values = PSC(3:end,jj)';
    
    % add the kernel field
    thePacket.kernel = kernelStruct;
    
    % add a metaData field
    thePacket.metaData = [];
    
    %% Perform the fit
    [paramsFit,fVal,modelResponseStruct] = ...
        temporalFit.fitResponse(thePacket,...
        'defaultParamsInfo', defaultParamsInfo, ...
        'searchMethod','linearRegression');
    
%     % loop over voxels --> returns a "cleaned" time series
%     for vxl = 1:size(meanSignal,1)
%         responseStruct.values = PSC(vxl,:,jj);
%         thePacket.response = responseStruct;
%         
%         % TFE linear regression here
%         
%     end
        
    % take mean across voxels
    
    
   % use TFE to run our model on the V1 mean 
        
end

    
    
    
    
    
    
    
    
    
