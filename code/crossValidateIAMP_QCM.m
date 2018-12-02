function [rmseMeanIAMP rmseMeanQCM rmseSemIAMP rmseSemQCM] = crossValidateIAMP_QCM(analysisParams,betas,timeCourseValues)
% Calculates the cross validated RMSE and SEM for both the QCM and IAMP
%
% Syntax:
%
%
% Description:
%
%
% Inputs:
%
%
% Outputs:
%    rmseMeanIAMP
%    rmseMeanQCM
%    rmseSemIAMP
%    rmseSemQCM
%
% Optional key/value pairs:
%    none

% Pick runs to leave out
leaveOut =[];
for ii = 1:size(betas,3)
    leaveOut = [leaveOut; randperm(size(betas,2))];
end



for jj = 1:size(betas,2)
    
    % Leave out run  
    for kk = 1:size(leaveOut,1)
        sheet = betas(:,:,kk);
        sheet(:,leaveOut(kk,jj)) =[];
        betasCV(:,:,kk) = sheet;
        
        % actual time course to compute rmse
        heldOutTimeCourse(:,kk) = timeCourseValues(:,leaveOut(kk,jj),kk);
    end
    
    baselineCV = betasCV(end,:,:);
    betasCV = betasCV(1:end-1,:,:);
    
    tmpBetas = [];
    for kk = 1:size(betasCV,3)
        tmpBetas =[tmpBetas; betasCV(:,:,kk)];
    end
    
    meanBetasCV = [mean(tmpBetas,2); mean(baselineCV(:))];
    
    %% Fit IAMP crfs with QCM
    % Set parameters and construct a QCM object.
    temporalFitQCM = tfeQCM('verbosity','none','dimension',analysisParams.theDimension);
    
    %% Set up contrast values matched to resoponse order
    % Set up stim order info to creat LMS contrast by timepoint matrix
    stimulusStruct.values   = [generateStimCombinations(analysisParams.contrastCoding,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension),[0;0]];
    stimulusStruct.timebase = 1:length(stimulusStruct.values);
    
    %% Snag response values from IAMP fit.
    thePacket.response.values = meanBetasCV';
    thePacket.response.timebase = 1:length(thePacket.response.values);
    
    %% Construct a packet for the QCM to fit.
    thePacket.stimulus = stimulusStruct;
    thePacket.kernel = [];
    thePacket.metaData = [];
    
    %% Fit
    % allow QCM to fit the offset
    defaultParamsInfo.noOffset = false;
    [paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket,'defaultParamsInfo',defaultParamsInfo);
    
    
    
    clear betasCV baselineCV
end

