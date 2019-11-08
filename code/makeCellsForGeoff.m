function [stimCells] = makeCellsForGeoff(subjId)
% This funciton makes a cell array contaiming the stimulus design matrix
% for each run. 
%
% Syntax:
%   [stimCells] = makeCellsForGeoff(subjId);
%
% Description:
%    Takes in a subject ID and returns a cell array {1xnumber of runs} that contains
%    the stimulus design matric for that run. 
%    {session 1 directions; session 2 directions; baseline condition;
%    attentional events}
%
% Inputs:
%    subjID                     - Subject ID used to get subject specific paramters
%                                 Can be 'LZ23', 'KAS25', 'AP26','LZ23_replication',
%                                 'KAS25_replication', 'AP26_replication'.
%                                 (string)
%
% Outputs:
%    analysisParams             - Returns analysisParams with any updates
%
% Optional key/value pairs:
%    none

% MAB 11/07/19

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);

% set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

stimCells = {};

for sessionNum = 1:length(analysisParams.sessionFolderName)
    
    % Set up Paths
    sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
    funcTextFile   = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'hcp','functionalRuns.txt');
    functionalPath = fullfile(sessionDir, 'hcp_func', analysisParams.sessionFolderName{sessionNum});
    
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    
    % Set up files.
    functionalRuns    = textFile2cell(funcTextFile);
    trialOrderFiles   = textFile2cell(trialOrderFile);
    functionalRuns    = fullfile(functionalPath,functionalRuns);
    
    % Gets the path to a text file that contains the mat file names needed
    % to get the trail order information for each run.
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderFiles = textFile2cell(trialOrderFile);
    
    % Number of acquisitions
    analysisParams.numAcquisitions = length(functionalRuns);
    
    % Save vars name
    saveName = [analysisParams.subjID,'_',analysisParams.sessionDate{sessionNum},'_area_V', num2str(analysisParams.areaNum),'_ecc_' num2str(analysisParams.eccenRange(1)) ,'_to_' ,num2str(analysisParams.eccenRange(2)) ,'_hcp.mat'];
    savePath = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'cleanTimeCourse');
    saveFullFile = fullfile(savePath,saveName);
    
    
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
        
        % get attention event regressor
        responseStruct.timeStep = analysisParams.timeStep;
        [~, eventsRegressor] = getAttentionEventTimes(block, responseStruct, 'timebase', thePacket.stimulus.timebase);
        
        % identify the data param file
        dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
        
        % We are about to load the data param file. First silence the warning
        % for EnumerableClassNotFound. Save the warning state.
        warningState = warning();
        warning('off','MATLAB:class:EnumerableClassNotFound')
        
        % Load and process the data param file
        load(dataParamFile);
        expParams = getExpParams(dataParamFile,analysisParams.TR,'hrfOffset', false, 'stripInitialTRs', false);
        
        stimBlocks =  createRegressors(expParams,analysisParams.baselineCondNum,totalTime,deltaT);
        
        % Place the session stim blocks(21xt) with the stimulus design
        % matrix(42xt)
        if sessionNum == 1
            tmpStimMat = [stimBlocks(1:end-1,:);zeros(size(stimBlocks(1:end-1,:)));stimBlocks(end,:);eventsRegressor];
            
        else
            tmpStimMat = [zeros(size(stimBlocks(1:end-1,:)));stimBlocks(1:end-1,:);stimBlocks(end,:);eventsRegressor];
        end
        
        stimCells{analysisParams.numAcquisitions*(sessionNum-1) + jj} = tmpStimMat;
        
    end
end
