function [paramFitsBootstrap,nrVals] = fitbootstrapQCM(analysisParams,iampParams,numIterations)

iampOBJ = tfeIAMP('verbosity','none');

for aa = 1:size(iampParams,1)
    for bb = 1:size(iampParams,2)
        iampParamsMatrix(:,bb+(aa-1)*size(iampParams,2)) = iampParams{aa,bb}.paramMainMatrix;
    end
end

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
    [qcmCrfMeanOBJ,qcmCrfMeanParams] = fitDirectionModel(analysisParams, 'qcmFit', {directionCrfMeanPacket});
    paramFitsBootstrap(:,iter) = qcmCrfMeanOBJ.paramsToVec(qcmCrfMeanParams{1});
    
    % Do some plotting of these fits
    nrVals(iter,:)= plotNakaRushtonFromParams(qcmCrfMeanParams{1}.crfAmp ,qcmCrfMeanParams{1}.crfExponent,qcmCrfMeanParams{1}.crfSemi,...
        'analysisParams',analysisParams,'plotFunction',false,'savePlot',false);
    
end


