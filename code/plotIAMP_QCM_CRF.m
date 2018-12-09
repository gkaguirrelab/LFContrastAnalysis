function [params] = plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit)
% This function plots the IAMP CRF and the IAMP-QCM CRF.
%
% Syntax:
%   [] = plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit)
%
% Description:
%    This function plots the IAMP fits and IAMP-QCM predictions from runIAMP_QCM.m as contrast response
%    functions (one plot per modulation direction).
%
% Inputs:
%    analysisParams            - Analysis parameter stuct set in analyzeLFContrast (Struct)
%    meanIAMPBetas             - Mean beta weughts across runs per contrast level and direction (vector)
%    semIAMPBetas              - Standard error of the beta weigths in meanIAMPBetas (vector)
%    paramsQCMFit              - Parameter fits to the QCM model (struct)
%
% Outputs:
%    params                    - Parameters of the Naka-Rushton fit. Each
%                                row of the matrix returned are the paramters
%                                fit to the IAMP CRF for a single run.
%                                The columns are:
%                                   Rmax  = params(1)
%                                   sigma = params(2)
%                                   n     = params(3)
%                                response = Rmax*[contrast^n]/[contrast^n + sigma^n]
%
% Optional key/value pairs:
%    none

% MAB 09/09/18

% Generate prediction to stimuli based on QCM fit to arbitrary stim
contrastSpacing  = linspace(max(analysisParams.contrastCoding),0,analysisParams.numSamples);
%contrastSpacing  = linspace(max(analysisParams.contrastCoding),min(analysisParams.contrastCoding),analysisParams.numSamples);
QCMStim.values   = generateStimCombinations(contrastSpacing,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension);
QCMStim.timebase = linspace(1,max(length(meanIAMPBetas(1:end-1))),length(QCMStim.values));
temporalFitQCM   = tfeQCM('verbosity','none','dimension',analysisParams.theDimension);
QCMResponses     = computeResponse(temporalFitQCM,paramsQCMFit,QCMStim,[]);

% Subplot size
rws = ceil(sqrt(size(analysisParams.directionCoding,2)));
cols = rws;
indx = length(analysisParams.contrastCoding);
figure
for ii = 1:size(analysisParams.directionCoding,2)
    
    % Get the contrast spacing for each plot.
    maxConVal = analysisParams.maxContrastPerDir(ii);
    maxContrastSpacing = contrastSpacing.*maxConVal;
    
    aa = (length(QCMResponses.values))/size(analysisParams.directionCoding,2);
    if ii == 1
        betas = meanIAMPBetas(1:indx);
        error = semIAMPBetas(1:indx);
        qcmSmooth = QCMResponses.values(1:aa);
    else
        betas = meanIAMPBetas((ii-1)*indx+1:ii*indx);
        error = semIAMPBetas((ii-1)*indx+1:ii*indx);
        qcmSmooth = QCMResponses.values((ii-1)*aa+1:ii*aa);
    end
    
    xAxis = [maxConVal.*analysisParams.contrastCoding, 0];
    betas = [betas;paramsQCMFit.crfOffset];
    error = [error;semIAMPBetas(end)];
    
    %% Plot the stuff
    subplot(rws,cols,ii); hold on
    p1 = errorbar(xAxis,betas,error,'k');
    p2 = plot(maxContrastSpacing,qcmSmooth,'r');
    
    % Plot Naka-Rushton Function
    [params(ii,:),f] = FittfeQCMComputeNakaRushton(xAxis',betas-paramsQCMFit.crfOffset);
    nrResponses = tfeQCMComputeNakaRushton(maxContrastSpacing,params(ii,2), params(ii,3),params(ii,1), paramsQCMFit.crfOffset);
    p3 = plot(maxContrastSpacing,nrResponses,'b');
    
    
    if isequal(analysisParams.directionCoding(:,ii),[1;1;0])
        pos = find(ismember(analysisParams.directionCoding',[1,1,0],'rows'));
        title(sprintf('L+M: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    elseif isequal(analysisParams.directionCoding(:,ii),[1;-1;0])
        pos = find(ismember(analysisParams.directionCoding',[1,-1,0],'rows'));
        title(sprintf('L-M: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    elseif isequal(analysisParams.directionCoding(:,ii),[1;0;0])
        pos = find(ismember(analysisParams.directionCoding',[1,0,0],'rows'));
        title(sprintf('L Isolating: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    elseif isequal(analysisParams.directionCoding(:,ii),[0;1;0])
        pos = find(ismember(analysisParams.directionCoding',[0,1,0],'rows'));
        title(sprintf('M Isolating: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    else
        title(sprintf('%s Degrees: Max Contrast = %s',num2str(analysisParams.LMVectorAngles(ii)),num2str(analysisParams.maxContrastPerDir(ii))))
    end
    ylabel('Mean Beta Weight')
    xlabel('Contrast')
    ylim([-0.8 1.1]);

    legend([p1, p2 p3], 'IAMP Model', 'QCM Fit', 'Naka-Rushton Fit', 'Location', 'northwest')
end

end
