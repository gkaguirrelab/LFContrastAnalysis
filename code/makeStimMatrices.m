function [stimCells] = makeStimMatrices(subjId,varargin)
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
%    analysisParams             - Returns analysisParams with any updates`
%
% Optional key/value pairs:
%    reorderCells               - Permute the order of of the cells based
%                                 on a vecotr input indication the postion
%                                 that current cell should be
%    padTheEnds                 - Add TRs assigned to the baseline
%                                 condition equal to the length of the MR
%                                 scan. (Default is to make the stimulus
%                                 design matrix the length specified by the
%                                 exp. paremters).
%

% MAB 11/07/19

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('subjId',@isstr);
p.addParameter('reorderCells',[],@isvector);
p.addParameter('padTheEnds',false,@islogical);
p.parse(subjId,varargin{:});

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
        if size(eventsRegressor,1) > 1 
            eventsRegressor = sum(eventsRegressor,1);
        end
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
        
        % Baseline pad the end
        if p.Results.padTheEnds
            if analysisParams.numClipFramesEnd == 0
                stimCells{analysisParams.numAcquisitions*(sessionNum-1) + jj} = tmpStimMat;
                clear tmpStimMat
            else
                tmpStimMat = [tmpStimMat zeros(size(tmpStimMat,1),analysisParams.numClipFramesEnd)];
                tmpStimMat(end - 1,analysisParams.expLengthTR:analysisParams.expLengthTR+analysisParams.numClipFramesEnd) = 1;
                stimCells{analysisParams.numAcquisitions*(sessionNum-1) + jj} = tmpStimMat;
                clear tmpStimMat
            end
        else
            stimCells{analysisParams.numAcquisitions*(sessionNum-1) + jj} = tmpStimMat;
            clear tmpStimMat
        end
    end
end

% reorder the cells to match the order they were run in ICA Fix
if ~isempty(p.Results.reorderCells)
    if length(p.Results.reorderCells) ~= length(stimCells)
        warning('Reordering vector must match the number of cells')
    end
    stimCells =stimCells(p.Results.reorderCells);
end
