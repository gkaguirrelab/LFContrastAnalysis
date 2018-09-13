%% analyzeLFContrast
%
% This script calls function in order to analyze the data for the
% LFContrast experiment.

% Convenience variables
projectName  = 'LFContrastAnalysis';
flywheelName = 'LFContrast';
subjID       = 'sub-HEROGKA1';
session      = 'ses-ResearchAguirre';
sessionFolderName = 'HERO_GKA1_2018-07-28';
sessionDir = fullfile(getpref('LFContrastAnalysis','projectRootDir'),sessionFolderName);

showPlots = false;

%% Relevant Nifti names for analysis
% NOTE: MB: What is the best way for us to set these files for analysis.
% especially when we have multiple session to combibe or multiple subjects
% functional runs
functionalRuns = {'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-2_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-4_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-4_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-5_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-5_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'};

confoundFiles  = {'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-1_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-2_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-2_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-3_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-3_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-4_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-4_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastPA_run-5_bold_confounds.tsv', ...
    'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-5_bold_confounds.tsv'};

numAcquisitions = length(functionalRuns);
% brain mask of function run for the reference volume in ANTs step
refFileName  = 'sub-HEROGKA1_ses-ResearchAguirre_task-tfMRILFContrastAP_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz';
% output files of Neuropythy (retinotopy template)
retinoFiles = {'HERO_gka1_native.template_angle.nii.gz','HERO_gka1_native.template_areas.nii.gz','HERO_gka1_native.template_eccen.nii.gz',};

% warp file name (product of running fmriprep)
warpFileName = 'sub-HEROGKA1_T1w_target-MNI152NLin2009cAsym_warp.h5';

% Set up paths to nifti and .h5 files
retinoPath     = fullfile(sessionDir,'neuropythy');
functionalPath = fullfile(sessionDir, 'fmriprep', subjID, session, 'func');
warpFilePath   = fullfile(sessionDir, 'fmriprep', subjID, session, 'anat');


%
% NOTE: MB: This section needs to be turned into a stand alone function that
% takes in a refrence image, a warp file and a movable image and saves out
% the warped image outpur
%

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
%
% NOTE:MB: Make function that takesin a functional run and a mask and
% returns a mean time series or voxel time courses. Talk with david on
% frisya meeting about packet making and saving
%
%Load mask nifti
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
trialOrderDir = fullfile(getpref(projectName,'melaDataPath'),'/Experiments/OLApproach_TrialSequenceMR/MRContrastResponseFunction/DataFiles/HERO_gka1/2018-07-28/session_1');
trialOrderFiles = {'CRF_session_1_scan1.mat' ...
    'CRF_session_1_scan2.mat', ...
    'CRF_session_1_scan3.mat', ...
    'CRF_session_1_scan4.mat', ...
    'CRF_session_1_scan5.mat', ...
    'CRF_session_1_scan6.mat', ...
    'CRF_session_1_scan7.mat', ...
    'CRF_session_1_scan8.mat', ...
    'CRF_session_1_scan9.mat', ...
    'CRF_session_1_scan10.mat'};

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
    
    %mean center the regressors
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
    %
    %NOTE:MB: Make a function that cleans up the time series
    %
    % loop over voxels --> returns a "cleaned" time series
    for vxl = 1:size(PSC,1)
        % place time series from this voxel into the packet
        thePacket.response.values = PSC(vxl,:);
        
        % TFE linear regression here
        [paramsFit, ~, QCMResponses] = temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
        confoundBetas(:,vxl) = paramsFit.paramMainMatrix;
        cleanRunData(vxl,:) = thePacket.response.values - QCMResponses.values;
    end
    
    % Store the mean across voxel confound values for this acquisition
    confoundBetasByAcq(:,jj) = mean(confoundBetas,2);
    
    % make stimulus values
    % Stim coding: 80% = 1, 40% = 2, 20% = 3, 10% = 4, 5% = 5, 0% = 6;
    baselineCondNum = 6;
    stimulusStruct.values =  createRegressors(expParams,baselineCondNum,totalTime,deltaT);
    
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
    thePacket.response.timebase = stimulusStruct.timebase;
    thePacket.response.values = median(cleanRunData,1);
    
    % add the kernel field
    thePacket.kernel = kernelStruct;
    
    % add a metaData field
    thePacket.metaData = [];
    
    %% Perform the fit
    [paramsFit,fVal,IAMPResponses] = ...
        temporalFit.fitResponse(thePacket,...
        'defaultParamsInfo', defaultParamsInfo, ...
        'searchMethod','linearRegression');
    
    temporalFit.plot(thePacket.response,'Color',[1 0 0]);
    temporalFit.plot(IAMPResponses,'Color',[0 1 0],'NewWindow',false);
    
    paramsFitIAMP{jj} = paramsFit;
    betas(:,jj)= paramsFit.paramMainMatrix;
    packetPocket{jj} = thePacket;
    modelResponses{jj} = IAMPResponses;
end

% Calculate mean and SEM of the betas
meanBetas = mean(betas,2);
semBeta = std(betas,0,2)./sqrt(numAcquisitions);
xPos = [100,50,25,12.5,6.25];

%% Save IAMP fit to time series data
save tempIAMPOutput
load tempIAMPOutput
%
%NOTE:MAB: below here is the QCM fit and 
%

%% Fit IAMP crfs with QCM
%
% Break out coefficients by stimulus color direction.
% Note that it is all hard coded, which we will need to 
% fix up at some point.
LminusMbetas = meanBetas(1:5)+ abs(meanBetas(21)); 
LplusMbetas = meanBetas(6:10)+abs(meanBetas(21));
LIsoBetas = meanBetas(11:15)+abs(meanBetas(21));
MIsoBetas = meanBetas(16:20)+abs(meanBetas(21));

% Set parameters and construct a QCM object.
theDimension = 2;
generatePlots = true;
temporalFitQCM = tfeQCM('verbosity','none','dimension',theDimension);

%% Set up contrast values matched to resoponse order
%
% Set up stim order info to creat LMS contrast by timepoint matrix
contrastCoding = [1, .5, .25, .125, .0625];
directionCoding = [1,1,1,0;-1,1,0,1;0,0,0,0]; %this 1 = L-M 2 = L+M 3 = L 4 = M;
maxContrastPerDir = [0.06,0.40,0.10,0.10]; % max contrast in the same order as above

% Now construct stimulus contrast description.
%
% Tag 0 contrast onto end as well.
stimulusStruct.values   = [generateStimCombinations(contrastCoding,directionCoding,maxContrastPerDir,theDimension),[0;0]];
stimulusStruct.timebase = 1:length(stimulusStruct.values);

%% Snag response values from IAMP fit.
thePacket.response.values = meanBetas(1:21)';
thePacket.response.timebase = 1:length(thePacket.response.values);
if (generatePlots)
    temporalFitQCM.plot(thePacket.response,'Color',[1 0 0]);
end

%% Construct a packet for the QCM to fit.
thePacket.stimulus = stimulusStruct;
thePacket.kernel = [];
thePacket.metaData = [];

%% Fit
[paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket);
fprintf('Model parameter from fits:\n');
temporalFitQCM.paramPrint(paramsQCMFit)

%% Generate prediction to stimuli based on QCM fit to arbitrary stim
numSamples = 25;
contrastSpacing = linspace(1,0.0625,numSamples);
QCMStim.values = [generateStimCombinations(contrastSpacing,directionCoding,maxContrastPerDir,theDimension),[0;0]];
QCMStim.timebase = linspace(1,max(thePacket.response.timebase),length(QCMStim.values));

QCMResponses = computeResponse(temporalFitQCM,paramsQCMFit,QCMStim,[]);

if (generatePlots)
    temporalFitQCM.plot(QCMResponses,'Color',[0 1 0],'NewWindow',false);
end

% NOTE: MB: Make a plot contrast response function. Should this be an 
%       update to the qcmTFE.plot method of a spereate function in LF 
%       Contrast analysis
%
%% Plot
%

if (generatePlots)
    figure
    subplot(2,2,1); hold on 
    error1 = semBeta(1:5);
    errorbar(xPos,LminusMbetas,error1)
    plot(contrastSpacing*100,QCMResponses.values(1:25)-QCMResponses.values(101))
    title('L-M: Max Contrast = 6%')
    ylabel('Mean Beta Weight')
    xlabel('Percent of Max Contrast')
    ylim([-0.2 1]);
    
    subplot(2,2,2); hold on
    error2 = semBeta(6:10);
    errorbar(xPos,LplusMbetas,error2)
    plot(contrastSpacing*100,QCMResponses.values(26:50)-QCMResponses.values(101))
    title('L+M: Max Contrast = 40%')
    ylabel('Mean Beta Weight')
    xlabel('Percent of Max Contrast')
    ylim([-0.2 1]);
    
    subplot(2,2,3); hold on
    error3 = semBeta(11:15);
    errorbar(xPos,LIsoBetas,error3)
    plot(contrastSpacing*100,QCMResponses.values(51:75)-QCMResponses.values(101))
    title('L Isolating: Max Contrast = 10%')
    ylabel('Mean Beta Weight')
    xlabel('Percent of Max Contrast')
    ylim([-0.2 1]);
    
    subplot(2,2,4); hold on
    error4 = semBeta(16:20);
    errorbar(xPos,MIsoBetas,error4)
    plot(contrastSpacing*100,QCMResponses.values(76:100)-QCMResponses.values(101))
    title('M Isolating: Max Contrast = 10%')
    ylabel('Mean Beta Weight')
    xlabel('Percent of Max Contrast')
    legend('IAMP CRF','QCM CRF')
    ylim([-0.2 1]);
end

LminusMcontrast = contrastCoding.*maxContrastPerDir(1);
LplusMcontrast= contrastCoding.*maxContrastPerDir(2);
LIsocontrast= contrastCoding.*maxContrastPerDir(3);
MIsocontrast= contrastCoding.*maxContrastPerDir(4);

IAMPBetas = {LminusMbetas,LplusMbetas,LIsoBetas,MIsoBetas};
contrastLevels = {LminusMcontrast,LplusMcontrast,LIsocontrast,MIsocontrast};
directions = {[1,-1],[1,1],[1,0],[0,1]};
thresh = 0.25;
hdl = plotIsorespContour(paramsQCMFit,IAMPBetas,contrastLevels,directions,thresh,[],'r');
thresh = 0.5;
hdl = plotIsorespContour(paramsQCMFit,IAMPBetas,contrastLevels,directions,thresh,hdl,'g');
thresh = 0.75;
hdl = plotIsorespContour(paramsQCMFit,IAMPBetas,contrastLevels,directions,thresh,hdl,'b');
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);

%% Now try plotting each run again, but also with mean and QCM fit params
for jj = 1:numAcquisitions   
    % Regenerate IAMP predictions to time series
    IAMPResponses = temporalFit.computeResponse(paramsFitIAMP{jj},packetPocket{jj}.stimulus,packetPocket{jj}.kernel);
    
    % Plot them
    temporalFit.plot(packetPocket{jj}.response,'Color',[1 0 0]);
    temporalFit.plot(IAMPResponses,'Color',[0 1 0],'NewWindow',false);
    
    % Doctor up the parameters to use mean IAMP values and plot again
    paramsFitIAMPMean = paramsFitIAMP{jj};
    paramsFitIAMPMean.paramMainMatrix(1:21) = meanBetas(1:21);
    IAMPResponsesMean = temporalFit.computeResponse(paramsFitIAMPMean,packetPocket{jj}.stimulus,packetPocket{jj}.kernel);
    temporalFit.plot(IAMPResponsesMean,'Color',[0 0.5 1],'NewWindow',false);
   
    % Doctor up parameters to use the QCM fit to the mean IAMP
    paramsFitIAMPQCM = paramsFitIAMP{jj};
    paramsFitIAMPQCM.paramMainMatrix(1:21) = fitResponseStructQCM.values(1:21');
    IAMPResponsesQCM = temporalFit.computeResponse(paramsFitIAMPQCM,packetPocket{jj}.stimulus,packetPocket{jj}.kernel);
    temporalFit.plot(IAMPResponsesQCM,'Color',[0 0 0],'NewWindow',false);
end
