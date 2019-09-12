function [qcmCrfMeanParams,nrVals] = fitbootstrapQCM(analysisParams,iampParams,numIterations)

iampOBJ = tfeIAMP('verbosity','none');
iampParamsMatrix = cell2mat(cellfun(@(x) x.paramMainMatrix, iampParams(:), 'UniformOutput', false)');

for iter = 1:numIterations
    bootIampParams = iampParams;
    for ii = 1:analysisParams.numSessions
        
        
        start = 1+(ii-1)*analysisParams.numAcquisitions;
        stop  = ii * analysisParams.numAcquisitions;
        tmpIampVals = iampParamsMatrix(:,start:stop);
        
        for jj = 1:size(iampParamsMatrix,1)
           indx = randi(10,1,10);
           tmpMatrix(jj,:,ii)= tmpIampVals(jj,indx);
        end
        
        for kk = 1:analysisParams.numAcquisitions
            bootIampParams{ii,kk}.paramMainMatrix= tmpMatrix(:,kk,ii);
        end
        
    end

    
    for ii = 1:analysisParams.numAcquisitions
        [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(bootIampParams(:,ii),'baselineMethod','averageBaseline');
        %[concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
    end
    averageIampParams = iampOBJ.averageParams(concatParams);
    
    directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,averageIampParams);
    
    %% Fit the direction based models to the mean IAMP beta weights
    %
    
    % Fit the CRF -- { } is because this expects a cell
    %[nrCrfOBJ,nrCrfParams] = fitDirectionModel(analysisParams, 'nrFit', {directionCrfMeanPacket});
    
    % Fit the CRF with the QCM -- { } is because this expects a cell
    [qcmCrfMeanOBJ,qcmCrfMeanParams{iter}] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});
    
    % Do some plotting of these fits
    tmpQCMParams = qcmCrfMeanParams{iter};
    [nrVals{iter}] = plotNakaRushtonFromParams(tmpQCMParams{1}.crfAmp ,tmpQCMParams{1}.crfExponent,tmpQCMParams{1}.crfSemi,...
        'analysisParams',analysisParams,'plotFunction',false,'savePlot',false);
    
end


