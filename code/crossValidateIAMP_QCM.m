function [rmseMeanIAMP, rmseMeanQCM, rmseSemIAMP, rmseSemQCM] = crossValidateIAMP_QCM(analysisParams,betas,timeCourseValues,paramsFitIAMP, packetPocket, varargin)
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
%    showPlots

p = inputParser;
p.addParameter('showPlots',true, @islogical);
p.parse(varargin{:});

% Pick runs to leave out
leaveOut =[];
for ii = 1:size(betas,3)
    leaveOut = [leaveOut; randperm(size(betas,2))];
end


% IAMP object
temporalFitIAMP = tfeIAMP('verbosity','none');

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
    betaLength = size(betasCV,1);
    
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
    
    fullTCHeldOut = [];
    fullIAMPpreds = [];
    fullQCMpreds  = [];
    for pp = 1:size(betasCV,3)
        
        heldOutRunNum = leaveOut(pp,jj) + size(betas,2)*(pp-1);
        % Doctor up the parameters to use mean IAMP values and plot again
        paramsFitIAMPMean = paramsFitIAMP{heldOutRunNum};
        paramsFitIAMPMean.paramMainMatrix(1:end-1) = [meanBetasCV(1+((pp-1)*betaLength):pp*betaLength);  mean(baselineCV(:))];
        IAMPResponsesMean = temporalFitIAMP.computeResponse(paramsFitIAMPMean,packetPocket{heldOutRunNum}.stimulus,packetPocket{heldOutRunNum}.kernel);
        
        % Doctor up parameters to use the QCM fit to the mean IAMP
        paramsFitIAMPQCM = paramsFitIAMP{heldOutRunNum};
        paramsFitIAMPQCM.paramMainMatrix(1:end-1) = [fitResponseStructQCM.values(1+((pp-1)*betaLength):pp*betaLength),  mean(baselineCV(:))]' ;
        IAMPResponsesQCM = temporalFitIAMP.computeResponse(paramsFitIAMPQCM,packetPocket{heldOutRunNum}.stimulus,packetPocket{heldOutRunNum}.kernel);
        
        
        if(analysisParams.generateCrossValPlots)
            originalTC = IAMPResponsesQCM;
            originalTC.values = heldOutTimeCourse(:,pp);
            temporalFitIAMP.plot(IAMPResponsesQCM,'Color',[1 0 0]);
            temporalFitIAMP.plot(IAMPResponsesMean,'Color',[0 1 0],'NewWindow',false);
            temporalFitIAMP.plot(originalTC,'Color',[0 0 0],'NewWindow',false);
        end
        
        % Calculate RMSE
        fullTCHeldOut = [fullTCHeldOut; heldOutTimeCourse(:,pp)] ;
        fullIAMPpreds = [fullIAMPpreds; IAMPResponsesMean.values'];
        fullQCMpreds  = [fullQCMpreds; IAMPResponsesQCM.values'];
        
        if pp ==2
            rmseIAMP(jj) = sqrt(mean((fullTCHeldOut - fullIAMPpreds).^2));
            rmseQCM(jj)  = sqrt(mean((fullTCHeldOut - fullQCMpreds).^2));
        end
    end
    clear betasCV baselineCV
end

rmseMeanIAMP = mean(rmseIAMP);
rmseMeanQCM  = mean(rmseQCM);
rmseSemIAMP  = std(rmseIAMP)./sqrt(length(rmseIAMP));
rmseSemQCM   = std(rmseQCM)./sqrt(length(rmseQCM));


if p.Results.showPlots
    figure; hold on
    c = categorical({'IAMP','QCM'});
    bar(c,[rmseMeanIAMP, rmseMeanQCM]);
    errorbar(1:2,[rmseMeanIAMP, rmseMeanQCM],[rmseSemIAMP,rmseSemQCM] ,'k.')
    ylim([0 1]);
    set(gca,'FontSize',12)
    axis square
end
