function [fig1] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampsPoints);
% This function plots the IAMP CRF and the IAMP-QCM CRF.
%
% Syntax:
%   [] = plotCF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit)
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
%    figHdl                    - Figure handle
%
% Optional key/value pairs:
%    none

% MAB 09/09/18

% Subplot size
rws = ceil(sqrt(size(analysisParams.directionCoding,2)));
cols = rws;
modelIndx = analysisParams.numSamples;
iampIndx = length(analysisParams.contrastCoding);
fig1 = figure;
contrastSpacing = crfStimulus.values(end,:);
fields = fieldnames(crfPlotStruct);

for ii = 1:size(analysisParams.directionCoding,2)
    
    for jj = 1:length(fields)
        
        theModelResp = eval(['crfPlotStruct.', fields{jj}]);
        
        % Get the contrast spacing for each plot.
        maxConVal = analysisParams.maxContrastPerDir(ii);
        
        if ii == 1
            crfValues = theModelResp.values(1:modelIndx);
            xAxisModels = contrastSpacing(1:modelIndx);
            iampVals = iampsPoints.paramMainMatrix(1:iampIndx)';
        else
            crfValues = theModelResp.values((ii-1)*modelIndx+1:ii*modelIndx);
            xAxisModels = contrastSpacing((ii-1)*modelIndx+1:ii*modelIndx);
            iampVals = iampsPoints.paramMainMatrix((ii-1)*iampIndx+1:ii*iampIndx)';
        end
        
        xAxisIamp = maxConVal.*analysisParams.contrastCoding;
        %% Plot the stuff
        subplot(rws,cols,ii); hold on
        p1 = plot(xAxisModels,crfValues,'color',theModelResp.color);
        scatter(xAxisIamp,iampVals,'k');
        
        ylabel('Mean Beta Weight')
        xlabel('Contrast')
        
        ylim([-0.3 1.4]);

    end
end


legend(feids)
end


