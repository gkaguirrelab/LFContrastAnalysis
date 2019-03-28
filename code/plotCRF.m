function [figHdl] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampsPoints)
% This function plots the IAMP CRF and the IAMP-QCM CRF.
%
% Syntax:
%   [figHdl] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampsPoints);
%
% Description:
%    This function plots the IAMP fits and IAMP-QCM predictions from runIAMP_QCM.m as contrast response
%    functions (one plot per modulation direction).
%
% Inputs:
%    analysisParams            - Analysis parameter stuct set in analyzeLFContrast (Struct)
%    crfPlotStruct             - A struct containing each model you want
%                                plotted as a field. Each model must subfields 
%                                of values (the CRF model predictions) and color 
%                                (the color values of the line)
%    crfStimulus               - The CRF stimulus used to make the model
%                                predictions 
%    iampsPoints               - The mean IAMP beta weights 
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

% indexind for models
modelIndx = analysisParams.numSamples;
iampIndx = length(analysisParams.contrastCoding);

% get x axis values 
contrastSpacing = crfStimulus.values(end,:);
fields = fieldnames(crfPlotStruct);

figHdl = figure; 

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
        p(jj) = plot(xAxisModels,crfValues,'color',theModelResp.color);
        q1    = scatter(xAxisIamp,iampVals,'*k');
        
        % put info 
        ylabel('Mean Beta Weight')
        xlabel('Contrast')
        title(sprintf('LM stim = %s', num2str(analysisParams.LMVectorAngles(ii))));
        ylim([-0.3 1.4]);

    end
end

fields{end+1}  = 'IAMP Points';
legend([p, q1], fields, 'Location','NorthWest')
set(gcf, 'Position',  [0, 0, 1800, 1300])
end


