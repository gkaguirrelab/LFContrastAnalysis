function [fullCleanData, analysisParams, voxelIndex] = getTimeCourse(analysisParams)
% Takes in the analysis params struct and returns a voxel by timepoint by
% aquisistion matrix for all the runs found fro the specicied session(s)
%
% Syntax:
%   [fullCleanData, analysisParams, voxelIndex] = getTimeCourse(analysisParams)
%
% Description:
%    This function takes in a struct that is specified in analyzeLFContrast.m
%    and returns a voxel by timepoint by aquisition matrix for all the
%    aquisitions specified in the text files housed in mela_analysis that
%    describe the session(s). f mutliple sessions, they will be
%    concatenated in the 3rd dimension
%
% Inputs:
%    analysisParams    - Stuct contianing relevenat info to the session that
%                        is defined in analyzeLFContrast.m. (string)
%
% Outputs:
%    fullCleanData       - The voxel by timepoint by aquisition matrix
%    analysisParams      - the input analysis params updated with the
%                          number of aquistidiond found per session
%    voxelIndex          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    none

% MAB 09/09/18


% Set up files and paths

% Initialize for output of various sessions.
fullCleanData = [];

% Loop over sessions
for sessionNum = 1:length(analysisParams.sessionFolderName)
    
    % Set up Paths
    sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
    funcTextFile   = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'fmriprep','functionalRuns.txt');
    confTexFile    = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'fmriprep','confounds.txt');
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    anatomyPath    = fullfile(sessionDir,'anatomy');
    retinoPath     = fullfile(anatomyPath,'neuropythy');
    functionalPath = fullfile(sessionDir, 'fmriprep', analysisParams.sessionFolderName{sessionNum}, 'fmriprep',  analysisParams.subjID, analysisParams.session{sessionNum}, 'func');
    warpFilePath   = fullfile(sessionDir, 'fmriprep', analysisParams.sessionFolderName{sessionNum},'fmriprep', analysisParams.subjID, 'anat');
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'melaDataPath'), analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    
    % Set up files.
    functionalRuns    = textFile2cell(funcTextFile);
    confoundFiles     = textFile2cell(confTexFile);
    trialOrderFiles   = textFile2cell(trialOrderFile);
    fullFileConfounds = fullfile(functionalPath,confoundFiles);
    functionalRuns    = fullfile(functionalPath,functionalRuns);
    refFile           = fullfile(functionalPath,analysisParams.refFileName);
    warpFile          = fullfile(warpFilePath,analysisParams.warpFileName);
    
    % Number of acquisitions
    analysisParams.numAcquisitions = length(functionalRuns);
    
    % Save vars name
    saveName = [analysisParams.subjID,'_',analysisParams.sessionDate{sessionNum},'_area_V', num2str(analysisParams.areaNum),'_ecc_' num2str(analysisParams.eccenRange(1)) ,'_to_' ,num2str(analysisParams.eccenRange(2)) ,'.mat'];
    savePath = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'cleanTimeCourse');
    saveFullFile = fullfile(savePath,saveName);
    
    % Load existing cleaned data
    if exist(saveFullFile)
        disp('cleaned time series file found')
        load(saveFullFile)
    else
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
        voxelTimeSeries = voxelTimeSeries(:,analysisParams.numClipFramesStart+1:end-analysisParams.numClipFramesEnd,:);
        
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
            load(dataParamFile);
            
            % restore warning state
            warning(warningState);
            
            % make stimulus timebase
            totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
            deltaT = analysisParams.TR*1000;
            thePacket.stimulus.timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
            thePacket.response.timebase = thePacket.stimulus.timebase;
            
            % get confound regressors
            confoundRegressors = getConfoundRegressors(fullFileConfounds{jj});
            
            % get attention event regressor
            responseStruct.timeStep = analysisParams.timeStep;
            [~, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', thePacket.stimulus.timebase);
            
            %mean center the regressors
            confoundRegressors = confoundRegressors - nanmean(confoundRegressors);
            confoundRegressors = confoundRegressors ./ nanstd(confoundRegressors);
            
            if size(confoundRegressors,1) > length(thePacket.stimulus.timebase)
                confoundRegressors = confoundRegressors(analysisParams.numClipFramesStart+1:end-analysisParams.numClipFramesEnd,:);
            end
            
            % Set up packet
            thePacket.kernel = [];
            thePacket.metaData = [];
            thePacket.stimulus.values = [confoundRegressors'; eventsRegressor];
            
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
                [paramsFit, ~, iampResponses] = temporalFit.fitResponse(thePacket,...
                    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
                confoundBetas(:,vxl) = paramsFit.paramMainMatrix;
                % Linear detrending of the timecourse
                cleanRunData(vxl,:,jj) = detrend(thePacket.response.values - iampResponses.values);
            end
        end
        %% Save out the clean time series brick
        save(saveFullFile,'cleanRunData','voxelIndex');
        
    end
    
    % Concatenate clean data across sessions.
    fullCleanData = cat(3,fullCleanData,cleanRunData);
    
end



