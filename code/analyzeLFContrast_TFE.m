%% analyzeLFContrast
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

%% Convenience variables
projectName  = 'LFContrastAnalysis';
flywheelName = 'LFContrast';
subjID       = 'sub-HEROgka1';
session      = 'ses-0411181853PM';
sessionFolderName = 'HERO_gka1_2018-04-11';
sessionDir = fullfile(getpref('LFContrastAnalysis','projectRootDir'),sessionFolderName);

%% Relevant Nifti names for analysis

% functional runs
functionalRuns = {'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

confoundFiles  = {'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-1_bold_confounds.tsv', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-1_bold_confounds.tsv', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-2_bold_confounds.tsv', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastAP_run-2_bold_confounds.tsv', ...
                  'sub-HEROgka1_ses-0411181853PM_task-tfMRILFContrastPA_run-3_bold_confounds.tsv'};
              
numAcquisitions = length(functionalRuns);
% brain mask of function run for the reference volume in ANTs step
refFileName  = 'sub-HEROgka1_ses-0411181853PM_task-tfMRIFLASHAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';

% output files of Neuropythy (retinotopy template)
retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};

% warp file name (product of running fmriprep)
warpFileName = 'sub-HEROgka1_ses-0411181853PM_T1w_target-MNI152NLin2009cAsym_warp.h5';

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
eccenRange  = [1 7];
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
[voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRuns,maskVol,'threshold', 0.5);

% Clip first two data points from the time series, as the signal is not yet
% steady state. We need to do more to either prevent these volumes from
% being saved, or to automatically detect this condition.
voxelTimeSeries = voxelTimeSeries(:,1:end,:);

%% Get trial order info:
trialOrderDir = fullfile(getpref(projectName,'melaDataPath'),'/Experiments/OLApproach_TrialSequenceMR/MRContrastResponseFunction/DataFiles/HERO_gka1/2018-04-11/session_1');
trialOrderFiles = {'CRF_session_1_scan1.mat', ...
                   'CRF_session_1_scan2.mat', ...
                   'CRF_session_1_scan3.mat', ...
                   'CRF_session_1_scan4.mat', ...
                   'CRF_session_1_scan5.mat'};
               
% make full file path to counfound tsv files
fullFileConfounds = fullfile(functionalPath,confoundFiles);

%% Construct the model object 
temporalFit = tfeIAMP('verbosity','none');

% Define the TR
    TR = 0.800;

%% Create a cell of stimulusStruct (one struct per run)
for jj = 1:numAcquisitions
    
    % identify the data param file
    dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});

    % We are about to load the data param file. First silence the warning
    % for EnumberableClass. Save the warning state.
    warningState = warning();
    warning('off','MATLAB:class:EnumerableClassNotFound')

    % Load and process the data param file
    load(dataParamFile);    
    expParams = getExpParams(dataParamFile,TR,'hrfOffset', false, 'stripInitialTRs', false);
    
    % restore warning state
    warning(warningState);
    
    % make stimulus timebase
    totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
    deltaT = 800;
    stimulusStruct.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
    responseStruct.timebase = stimulusStruct.timebase;
    
    % get confound regressors 
    confoundRegressors = getConfoundRegressors(fullFileConfounds{jj});

    % trim the first two frames, and normalize the regressors
    confoundRegressors = confoundRegressors(1:end,:);
    confoundRegressors = confoundRegressors - nanmean(confoundRegressors);
    confoundRegressors = confoundRegressors ./ nanstd(confoundRegressors);
    
    thePacket.kernel = [];
    thePacket.metaData = [];
    thePacket.stimulus.timebase = stimulusStruct.timebase;
    thePacket.stimulus.values = confoundRegressors';
    
    defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);
    % get the data for all masked voxel in a run 
    runData = voxelTimeSeries(:,:,jj);
    
    % convert to percent signal change relative to the mean
    voxelMeanVec = mean(runData,2);
    PSC = 100*((runData - voxelMeanVec)./voxelMeanVec);

    % timebase will be the same for every voxel
    thePacket.response.timebase = stimulusStruct.timebase;

    % loop over voxels --> returns a "cleaned" time series
    for vxl = 1:size(PSC,1)
        % place time series from this voxel into the packet
        thePacket.response.values = PSC(vxl,:);
        
        % TFE linear regression here
        [paramsFit, ~, modelResponseStruct] = temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
        confoundBetas(:,vxl) = paramsFit.paramMainMatrix;
        cleanRunData(vxl,:) = thePacket.response.values - modelResponseStruct.values;
    end

    % Store the mean across voxel confound values for this acquisition
    confoundBetasByAcq(:,jj) = mean(confoundBetas,2);
    
    % make stimulus values
    % Stim coding: 80% = 1, 40% = 2, 20% = 3, 10% = 4, 5% = 5, 0% = 6;
    stimulusStruct.values = zeros(length(unique(expParams(:,3))),totalTime/deltaT);
    for kk = 1:size(expParams,1)
        stimulusStruct.values(expParams(kk,3),expParams(kk,1):expParams(kk,2)) = 1;
    end
    
    %[ * NOTE: MB: make sure the timestep is loaded from the pulse params
    %istead of set here]
    responseStruct.timeStep = 1/100;
    
    % get attention event regressor 
    [attentionEventTimes, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', thePacket.response.timebase);
    
    % add attention events to regressor matrix 
    stimulusStruct.values = [stimulusStruct.values;eventsRegressor];
    
    % Set the number of instances.
    defaultParamsInfo.nInstances = size(stimulusStruct.values,1);

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
    thePacket.stimulus.timebase = stimulusStruct.timebase;
    thePacket.stimulus.values   = stimulusStruct.values;
    
    % add the response field
    thePacket.response.timebase =stimulusStruct.timebase;
    thePacket.response.values = median(cleanRunData,1);
    
    % add the kernel field
    thePacket.kernel = kernelStruct;
    
    % add a metaData field
    thePacket.metaData = [];
    
    %% Perform the fit
    [paramsFit,fVal,modelResponseStruct] = ...
        temporalFit.fitResponse(thePacket,...
        'defaultParamsInfo', defaultParamsInfo, ...
        'searchMethod','linearRegression');
    
    betas(:,jj)= paramsFit.paramMainMatrix;
    packetPocket{jj} = thePacket;
    modelResponses{jj} = modelResponseStruct;
end

%% Plot the CRF

meanBetas = mean(betas,2);
semBeta = std(betas,0,2)./sqrt(numAcquisitions);
contrastVals = [80, 40, 20,10,5,0];
figure; hold on
errorbar(contrastVals,meanBetas(1:6),semBeta(1:6))
title('Contrast Response Function')
xlabel('Contrast Level')
ylabel('Mean Beta Weight')
    
    
    
    
    
    
    
