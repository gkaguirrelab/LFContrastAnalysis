function [params,fullCleanData] = simulateDataFromExpParams(analysisParams,betaWeights,numDirections,numContrast,numVoxels,varargin)
% Generates simulated voxel timecourses based on betaweights and 
%
% Syntax:
%    thePackets = generateSamplePackets(betaWeights,numDirections,numContrast,numPackets)
%
% Description:
%    This function takes in a vector of beta weights (one for each
%    regressor) and returns a random time course repsonse packet.
%
% Inputs:
%    analysisParams            - Beta weight values for each condtion
%    betaWeights               -  Number of directions
%    numDirections             - Number of contrast conditions
%    numContrast               - Number of packets to be generated
%    numVoxels                 -
% Outputs:
%    thePackets                - Packets of simulated time courses
% Optional key/value pairs:
%    baseline                  - Baseline added to each run
%    linDetrending             - 

% MAB 06/10/19

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isvector);
p.addRequired('betaWeights',@isnumeric);
p.addRequired('numDirections',@isnumeric);
p.addRequired('numContrast',@isnumeric);
p.addRequired('numVoxels',@isnumeric);
p.addParameter('baseline',0,@isnumeric);
p.addParameter('linDetrending',false,@islogical)

p.parse(analysisParams, betaWeights,numDirections,numContrast,numVoxels,varargin{:});

counter = 1;
for sessionNum = 1:analysisParams.numSessions
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderFiles = textFile2cell(trialOrderFile);
    
    
    for jj = 1:analysisParams.numAcquisitions
        dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
        expParams = getExpParams(dataParamFile,analysisParams.TR,'hrfOffset', false, 'stripInitialTRs', false);
        
        [params{counter}, data] = generateSampleVoxels(betaWeights,numDirections,numContrast,numVoxels, 'realExpParams',expParams);
        
        fullCleanData(:,:,counter) = data + p.Results.baseline;
        
        if p.Results.linDetrending
            fullCleanData(:,:,counter) = detrend(data')';
        end
        
        counter = counter+1;
    end
end

