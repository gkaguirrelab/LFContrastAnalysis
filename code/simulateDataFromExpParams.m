function [params,fullCleanData] = simulateDataFromExpParams(analysisParams,betaWeights,numDirections,numContrast,numVoxels)

counter = 1;
for sessionNum = 1:analysisParams.numSessions
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderFiles = textFile2cell(trialOrderFile);
    
    
    for jj = 1:analysisParams.numAcquisitions
        dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
        expParams = getExpParams(dataParamFile,analysisParams.TR,'hrfOffset', false, 'stripInitialTRs', false);
        
        [params{counter}, data] = generateSampleVoxels(betaWeights,numDirections,numContrast,numVoxels, 'realExpParams',expParams);
        
        %fullCleanData(:,:,counter) = detrend(data')';
        fullCleanData(:,:,counter) = data - 0.5;
        
        counter = counter+1;
    end
end

