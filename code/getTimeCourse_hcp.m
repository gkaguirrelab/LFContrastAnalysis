function [fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams,varargin)
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
%    regressAttenEvents  - Option to either regress out attentional events
%                          or keep them in.
%    polyFitOrder        - Order of the polynomial fit for trend removal.
% MAB 09/09/18

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addParameter('regressAttenEvents',true,@islogical);
p.addParameter('polyFitOrder',5,@isnumeric);
p.addParameter('highpass',false,@islogical);
p.addParameter('wholeBrain',false,@islogical);
p.parse(analysisParams,varargin{:});

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
    
    savePathROI  = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis','MNI_ROIs');
    
    if analysisParams.useSubcortROI
        roiSaveName        = analysisParams.subcortROI;
    elseif p.Results.wholeBrain
        voxelsSaveName     = 'wholeBrain';
        roiSaveName        = [voxelsSaveName,'.dscalar.nii'];
    else
        voxelsSaveName     = ['V', num2str(analysisParams.areaNum), '_', analysisParams.hemisphere, '_ecc_', num2str(analysisParams.eccenRange(1)), '_to_', num2str(analysisParams.eccenRange(2))];
        roiSaveName        = [voxelsSaveName,'.dscalar.nii'];
    end
    maskFullFile = fullfile(savePathROI,roiSaveName);
    
    % Number of acquisitions
    analysisParams.numAcquisitions = length(functionalRuns);
    
    % Save vars name
    if analysisParams.useSubcortROI
        [~,tmp,~] = fileparts(analysisParams.subcortROI);
        [~,regionName,~] = fileparts(tmp);
        dataSaveName =[analysisParams.subjID,'_',analysisParams.sessionDate{sessionNum},'_area_',regionName, '_hcp.mat'];
    elseif p.Results.wholeBrain
        dataSaveName = [analysisParams.subjID,'_',analysisParams.sessionDate{sessionNum},'_wholeBrain_hcp.mat'];
    else
        dataSaveName = [analysisParams.subjID,'_',analysisParams.sessionDate{sessionNum},'_area_V', num2str(analysisParams.areaNum),'_ecc_' num2str(analysisParams.eccenRange(1)) ,'_to_' ,num2str(analysisParams.eccenRange(2)) ,'_hcp.mat'];
    end
    savePath = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'cleanTimeCourse');
    saveFullFile = fullfile(savePath,dataSaveName);
    
    % Load existing cleaned data
    saveFileStatus = 0;
    if exist(saveFullFile) && exist(maskFullFile)
        saveFileStatus = 1;
        disp('cleaned time series file found')
        load(saveFullFile)
        censorPoints{sessionNum,:} = saveCPoints;
        [ maskMatrix ] = loadCIFTI(maskFullFile,'workbenchPath',getpref(analysisParams.projectName,'wbPath'));
        voxelIndex = find(maskMatrix);
        
    else
        
        %% Create restricted V1 mask
        if p.Results.wholeBrain
            maskMatrix = ones(91282,1);
            voxelIndex = find(maskMatrix);
        elseif exist(maskFullFile)
            display(sprintf('ROI Found: %s',roiSaveName))
            [ maskMatrix ] = loadCIFTI(maskFullFile,'workbenchPath',getpref(analysisParams.projectName,'wbPath'));
            voxelIndex = find(maskMatrix);
        else
            %melaAnalysisPath = '/Users/michael/labDropbox/MELA_analysis/';
            pathToBensonMasks = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'), 'mriTOMEAnalysis','flywheelOutput','benson');
            pathToBensonMappingFile = fullfile(pathToBensonMasks,'indexMapping.mat');
            pathToTemplateFile = fullfile(pathToBensonMasks,'template.dscalar.nii');
            
            maskMatrix = makeMaskFromRetinoCIFTI(analysisParams.areaNum, analysisParams.eccenRange, analysisParams.anglesRange, analysisParams.hemisphere, ...
                'saveName', maskFullFile, 'pathToBensonMasks', pathToBensonMasks, 'pathToTemplateFile', pathToTemplateFile, ...
                'pathToBensonMappingFile', pathToBensonMappingFile, 'threshold', analysisParams.threshold);
            voxelIndex = find(maskMatrix);
        end
        
        %% Extract Signal from voxels
        if analysisParams.useSubcortROI
            voxelsSaveName = regionName;
        elseif p.Results.wholeBrain
            voxelsSaveName = wholeBrain;
        end
        saveVoxelTimeSeriesName = fullfile(functionalPath,'tfMRI_LFContrast_AllRuns',['voxelTimeSeries_' voxelsSaveName '.mat']);
        if exist(saveVoxelTimeSeriesName)
            theVars = load(saveVoxelTimeSeriesName);
            voxelTimeSeries = theVars.voxelTimeSeries;
        else
            for ii = 1:length(functionalRuns)
                % load nifti for functional run
                cifti = loadCIFTI(functionalRuns{ii},'workbenchPath',getpref(analysisParams.projectName,'wbPath'));
                voxelTimeSeries(:,:,ii) = cifti(logical(maskMatrix),:);
                
            end
            save(saveVoxelTimeSeriesName,'voxelTimeSeries')
        end
        
        % Clip initial frames if specified
        numTimePoints = size(voxelTimeSeries,2);
        
        
        if numTimePoints > analysisParams.expLengthTR
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
            
            % Get the censored time-point times
            [cPoints{sessionNum,jj}, percentCensored] = findCensoredPoints(analysisParams,movementRegressorsFull(:,1:6),...
                'plotMotion',false, 'distMetric', 'l2','addBuffer',[1,1]);
            
            relativeMovementRegressors = movementRegressorsFull(:,7:12);
            
            % get attention event regressor
            responseStruct.timeStep = analysisParams.timeStep;
            [~, eventDeltas] = getAttentionEventTimes(block, responseStruct, 'timebase', thePacket.stimulus.timebase);
            
            % Convolve with hrf
            eventStruct.values = eventDeltas;
            eventStruct.timebase = timeBase;
            outputStruct = applyKernel(temporalFit,eventStruct,analysisParams.HRF);
            eventsRegressor = outputStruct.values;
            
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
            if p.Results.regressAttenEvents
                if unique(eventsRegressor) == 0
                    thePacket.stimulus.values = [confoundRegressors];
                else
                    thePacket.stimulus.values = [confoundRegressors; eventsRegressor];
                end
            else
                thePacket.stimulus.values = [confoundRegressors];
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
                f = fit(timeBase',S',['poly' num2str(p.Results.polyFitOrder)], 'Exclude', find(isnan(cTimeSeries)));
                cleanRunData(vxl,:,jj) = S - f(timeBase)';
                
            end
            clear thePacket
        end
        %% Save out the clean time series brick
        saveCPoints = cPoints(sessionNum,:);
        save(saveFullFile,'cleanRunData','maskMatrix','saveCPoints');
        clear voxelTimeSeries
    end
    
    % Concatenate clean data across sessions.
    fullCleanData = cat(3,fullCleanData,cleanRunData);
    
end

if saveFileStatus == 0
    analysisParams.censorPoints = cPoints;
else
    analysisParams.censorPoints = censorPoints{1,:};
    for ll = 2:length(analysisParams.sessionFolderName)
        analysisParams.censorPoints = [analysisParams.censorPoints;censorPoints{ll,:}];
    end
end

