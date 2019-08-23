function [fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams)
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
    confTexFile    = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'hcp','confounds.txt');
    funcTextFile   = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'hcp','functionalRuns.txt');
    functionalPath = fullfile(sessionDir, 'hcp_func', analysisParams.sessionFolderName{sessionNum});
    
    
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    anatomyPath    = fullfile(sessionDir,'anatomy');
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    
    % Set up files.
    functionalRuns    = textFile2cell(funcTextFile);
    confoundFiles     = textFile2cell(confTexFile);
    trialOrderFiles   = textFile2cell(trialOrderFile);
    functionalRuns    = fullfile(functionalPath,functionalRuns);
    fullFileConfounds = fullfile(functionalPath,confoundFiles);
    
    % Number of acquisitions
    analysisParams.numAcquisitions = length(functionalRuns);
    
    % Save vars name
    saveName = [analysisParams.subjID,'_',analysisParams.sessionDate{sessionNum},'_area_V', num2str(analysisParams.areaNum),'_ecc_' num2str(analysisParams.eccenRange(1)) ,'_to_' ,num2str(analysisParams.eccenRange(2)) ,'_hcp.mat'];
    savePath = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'cleanTimeCourse');
    saveFullFile = fullfile(savePath,saveName);
    
    % Load existing cleaned data
    if exist(saveFullFile)
        disp('cleaned time series file found')
        load(saveFullFile)
    else
        
        %% Create restricted V1 mask
        savePathROI  = fullfile(getpref(analysisParams.projectName,'projectRootDir'),'MNI','ROIs');
        saveName     = ['V', num2str(analysisParams.areaNum), '_', analysisParams.hemisphere, '_ecc_', num2str(analysisParams.eccenRange(1)), '_to_', num2str(analysisParams.eccenRange(2)),'.dscalar.nii'];
        maskFullFile = fullfile(savePathROI,saveName);
        
        if exist(maskFullFile)
            [ maskMatrix ] = loadCIFTI(maskFullFile);
        else
            %melaAnalysisPath = '/Users/michael/labDropbox/MELA_analysis/';
            pathToBensonMasks = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');
            pathToBensonMappingFile = fullfile(pathToBensonMasks,'indexMapping.mat');
            pathToTemplateFile = fullfile(pathToBensonMasks,'template.dscalar.nii');
            
            maskMatrix = makeMaskFromRetinoCIFTI(analysisParams.areaNum, analysisParams.eccenRange, analysisParams.anglesRange, analysisParams.hemisphere, ...
                'saveName', maskFullFile, 'pathToBensonMasks', pathToBensonMasks, 'pathToTemplateFile', pathToTemplateFile, ...
                'pathToBensonMappingFile', pathToBensonMappingFile, 'threshold', analysisParams.threshold);
        end
        
        %% Extract Signal from voxels
        saveVoxelTimeSeriesName = fullfile(functionalPath,'tfMRI_LFContrast_AllRuns','voxelTimeSeries.mat');
        if exist(saveVoxelTimeSeriesName)
            theVars = load(saveVoxelTimeSeriesName);
            voxelTimeSeries = theVars.voxelTimeSeries;
        else
            for ii = 1:length(functionalRuns)
                % load nifti for functional run
                cifti = loadCIFTI(functionalRuns{ii});
                voxelTimeSeries(:,:,ii) = cifti(logical(maskMatrix),:);
                
            end
            save(saveVoxelTimeSeriesName,'voxelTimeSeries')
        end
        
        % Clip initial frames if specified
        numTimePoints = size(voxelTimeSeries,2);
        
        
        if numTimePoints > 360
            voxelTimeSeries = voxelTimeSeries(:,analysisParams.numClipFramesStart+1:end-analysisParams.numClipFramesEnd,:);
        end
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
            timeBase = linspace(0,totalTime-deltaT,totalTime/deltaT);
            
            thePacket.stimulus.timebase = timeBase;
            thePacket.response.timebase = timeBase;
            
            % get confound regressors
            %Movement_Regressors.txt
            fields_per_line = 12;
            fileID = fopen(fullFileConfounds{jj},'r');
            formatSpec = '%f';
            textVector = fscanf(fileID,formatSpec);
            fclose(fileID);
            movementRegressorsFull     = reshape(textVector,[fields_per_line,numTimePoints])';
            [cPoints{sessionNum,jj}, percentCensored] = findCensoredPoints(movementRegressorsFull,'plotMotion',false, 'distMetric', 'l2');
            relativeMovementRegressors = movementRegressorsFull(:,7:12);
            
            % get attention event regressor
            responseStruct.timeStep = analysisParams.timeStep;
            [~, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', thePacket.stimulus.timebase);
            
            %mean center the regressors
            relativeMovementRegressors = relativeMovementRegressors - nanmean(relativeMovementRegressors);
            relativeMovementRegressors = relativeMovementRegressors ./ nanstd(relativeMovementRegressors);
            
            %check for nans
            nanCol = find(all(isnan(relativeMovementRegressors),1))
            
            if ~isempty(nanCol)
                relativeMovementRegressors(:,nanCol) = [];
            end
            
            
            if size(relativeMovementRegressors,1) > length(thePacket.stimulus.timebase)
                confoundRegressors = relativeMovementRegressors(analysisParams.numClipFramesStart+1:end-analysisParams.numClipFramesEnd,:)';
            else
                confoundRegressors =relativeMovementRegressors';
            end
            
            % Set up packet
            thePacket.kernel = [];
            thePacket.metaData = [];
            if unique(eventsRegressor) == 0
                thePacket.stimulus.values = [confoundRegressors];
            else
                thePacket.stimulus.values = [confoundRegressors; eventsRegressor];
            end
            
            defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);
            % get the data for all masked voxel in a run
            runData = voxelTimeSeries(:,:,jj);
            
            % convert to percent signal change relative to the mean
            voxelMeanVec = mean(runData,2);
            PSC = 100*((runData - voxelMeanVec)./voxelMeanVec);
            
            sprintf('session %s, run %s', num2str(sessionNum), num2str(jj))
            
            % loop over voxels --> returns a "cleaned" time series
            thePacket.stimulus.values(:,cPoints{sessionNum,jj}) = [];
            thePacket.stimulus.timebase(cPoints{sessionNum,jj}) = [];
            thePacket.response.timebase(cPoints{sessionNum,jj}) = [];
            for vxl = 1:size(PSC,1)
                % place time series from this voxel into the packet
                thePacket.response.values = PSC(vxl,:);
                cTimeSeries =  censorFrames(cPoints{sessionNum,jj},thePacket.response.values);
                
                % Remove censor Points for regression
                thePacket.response.values(cPoints{sessionNum,jj}) =[];

                % TFE linear regression here
                if any(isnan(thePacket.response.values))
                    [paramsFit, ~, iampResponses] = temporalFit.fitResponse(thePacket,...
                        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','fmincon');
                else
                    [paramsFit, ~, iampResponses] = temporalFit.fitResponse(thePacket,...
                        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
                end
                % add back NaNs where censored points were removed
                nanVec = nan(size(timeBase));
                tmp = ones(size(timeBase));
                tmp(cPoints{sessionNum,jj}) = 0;
                nanVec(tmp==1) = iampResponses.values;
                
                % Linear detrending of the timecourse
                S = cTimeSeries - nanVec;
                f = fit(timeBase',S','poly1', 'Exclude', find(isnan(cTimeSeries)));
                cleanRunData(vxl,:,jj) = S - f(timeBase)';
            end
            clear thePacket
        end
        %% Save out the clean time series brick
        save(saveFullFile,'cleanRunData','maskMatrix');
        clear voxelTimeSeries
    end
    
    % Concatenate clean data across sessions.
    fullCleanData = cat(3,fullCleanData,cleanRunData);
    
end



