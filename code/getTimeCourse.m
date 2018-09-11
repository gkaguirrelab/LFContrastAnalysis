function [cleanRunData, analysisParams, voxelIndex] = getTimeCourse(analysisParams)
% Takes in a text file name and retuns a cell of the lines of the text file
%
% Syntax:
%   filesCell = textFile2cell(inFile)
%
% Description:
%    This function takes in a file name for a trext file and returns a cell
%    that is composed of the lines of the text file. Example of this would
%    be a text file of file names the output is a cell of files names.
%
% Inputs:
%    inFile            - File name of a text file. (string)
%
% Outputs:
%    fileCell          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    none

% MAB 09/09/18


% set up files and paths
sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.sessionFolderName);
funcTextFile   = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName,'fmriprep','functionalRuns.txt');
confTexFile    = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName,'fmriprep','confounds.txt');
trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName,'experimentFiles','dataFiles.txt');
retinoPath     = fullfile(sessionDir,'neuropythy');
functionalPath = fullfile(sessionDir, 'fmriprep', analysisParams.subjID, analysisParams.session, 'func');
warpFilePath   = fullfile(sessionDir, 'fmriprep', analysisParams.subjID, analysisParams.session, 'anat');
trialOrderDir  = fullfile(getpref(analysisParams.projectName,'melaDataPath'), analysisParams.expSubjID,analysisParams.sessionDate,analysisParams.sessionNumber);

functionalRuns  = textFile2cell(funcTextFile);
confoundFiles   = textFile2cell(confTexFile);
trialOrderFiles = textFile2cell(trialOrderFile);
analysisParams.numAcquisitions = length(functionalRuns);

fullFileConfounds = fullfile(functionalPath,confoundFiles);
functionalRuns    = fullfile(functionalPath,functionalRuns);
refFile           = fullfile(functionalPath,analysisParams.refFileName);
warpFile          = fullfile(warpFilePath,analysisParams.warpFileName);

%% Create restricted V1 mask
% load ecc nifti file
eccenPos       = find(~cellfun(@isempty,strfind(analysisParams.retinoFiles,'eccen')));
[~,tempName,~] = fileparts(analysisParams.retinoFiles{eccenPos});
[~,outName,~]  = fileparts(tempName);
eccenFileName  = fullfile(retinoPath,[outName '.nii.gz']);
eccen          = MRIread(eccenFileName);

% load areas nifti file
areasPos       = find(~cellfun(@isempty,strfind(analysisParams.retinoFiles,'areas')));
[~,tempName,~] = fileparts(analysisParams.retinoFiles{areasPos});
[~,outName,~]  = fileparts(tempName);
areasFileName  = fullfile(retinoPath,[outName,'.nii.gz']);
areas          = MRIread(areasFileName);

% make mask from the area and eccentricity maps
[~,maskSaveName] = makeMaskFromRetino(eccen,areas,analysisParams.areaNum,analysisParams.eccenRange,retinoPath);
files2warp = {'HERO_gka1_T1.nii.gz',maskSaveName};

%% Apply the warp to the mask and T1 files using ANTs
inFiles = fullfile(retinoPath,files2warp);
applyANTsWarpToData(inFiles, warpFile, refFile);

%% Extract Signal from voxels
maskPos         = find(~cellfun(@isempty,strfind(files2warp,'mask')));
[~,tempName,~]  = fileparts(files2warp{maskPos});
[~,tmpName,~]   = fileparts(maskSaveName);
[~,outName,~]   = fileparts(tmpName);
maskOutFileName = fullfile(retinoPath,[outName '_MNI_resampled.nii.gz']);
mask            = MRIread(maskOutFileName);
maskVol         = mask.vol;

% extract the mean signal from voxels
[voxelTimeSeries, voxelIndex] = extractTimeSeriesFromMask(functionalRuns,maskVol,'threshold', 0.5);

% Clip initial frames if specified
voxelTimeSeries = voxelTimeSeries(:,analysisParams.numClipFrames+1:end,:);

%% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

%% Create a cell of stimulusStruct (one struct per run)
for jj = 1:analysisParams.numAcquisitions
    
    % identify the data param file
    dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
    
    % We are about to load the data param file. First silence the warning
    % for EnumberableClass. Save the warning state.
    warningState = warning();
    warning('off','MATLAB:class:EnumerableClassNotFound')
    
    % Load and process the data param file
    load(dataParamFile,'protocolParams');
    
    % restore warning state
    warning(warningState);
    
    % make stimulus timebase
    totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
    deltaT = analysisParams.TR*1000;
    thePacket.stimulus.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
    thePacket.response.timebase = thePacket.stimulus.timebase;
    
    % get confound regressors
    confoundRegressors = getConfoundRegressors(fullFileConfounds{jj});
    
    %mean center the regressors
    confoundRegressors = confoundRegressors - nanmean(confoundRegressors);
    confoundRegressors = confoundRegressors ./ nanstd(confoundRegressors);
    
    thePacket.kernel = [];
    thePacket.metaData = [];
    thePacket.stimulus.values = confoundRegressors';
    
    defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);
    % get the data for all masked voxel in a run
    runData = voxelTimeSeries(:,:,jj);
    
    % convert to percent signal change relative to the mean
    voxelMeanVec = mean(runData,2);
    PSC = 100*((runData - voxelMeanVec)./voxelMeanVec);
    
    % loop over voxels --> returns a "cleaned" time series
    for vxl = 1:size(PSC,1)
        % place time series from this voxel into the packet
        thePacket.response.values = PSC(vxl,:);
        
        % TFE linear regression here
        [paramsFit, ~, QCMResponses] = temporalFit.fitResponse(thePacket,...
            'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
        confoundBetas(:,vxl) = paramsFit.paramMainMatrix;
        cleanRunData(vxl,:,jj) = thePacket.response.values - QCMResponses.values;
    end  
end

%% Save out the clean time series brick
saveName = [analysisParams.subjID,'_',analysisParams.sessionDate,'_area_V', num2str(analysisParams.areaNum),'_ecc_' num2str(analysisParams.eccenRange(1)) ,'_to_' ,num2str(analysisParams.eccenRange(2)) ,'.mat'];
savePath = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName,'cleanTimeCourse');
saveFullFile = fullfile(savePath,saveName);
save(saveFullFile,'cleanRunData','voxelIndex','cleanRunData');

