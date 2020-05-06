function [figHdl] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampPoints, iampError, varargin)
% This function plots the IAMP CRF and the IAMP-QCM CRF.
%
% Syntax:
%   [figHdl] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampsPoints);
%
% Description:
%    This function plots the IAMP fits and IAMP-QCM predictions from
%    runIAMP_QCM.m as contrast response functions (one plot per modulation
%    direction).
%
% Inputs:
%    analysisParams            - Analysis parameter stuct set in
%                                analyzeLFContrast (Struct)
%    crfPlotStruct             - A struct containing each model you want
%                                plotted as a field. Each model must
%                                subfields of values (the CRF model
%                                predictions) and color (the color values
%                                of the line)
%    crfStimulus               - The CRF stimulus used to make the model
%                                predictions
%    iampsPoints               - The mean IAMP beta weights
%    iampError                 - Error bars for iampPoints
%
% Outputs:
%    figHdl                    - Figure handle
%
% Optional key/value pairs:
%    subtractBaseline          - The baseline value to be subtracted from
%                                the CRF values
%    iampColor                 - Optional color for IAMP markers
%    indivBootCRF              - Plot each draw of a bootstrap analysis in
%                                each corresponding CRF

% MAB 09/09/18

% Subplot size
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addRequired('crfPlotStruct',@isstruct);
p.addRequired('crfStimulus',@isstruct);
p.addRequired('iampPoints',@isstruct);
p.addRequired('iampError',@isstruct);
p.addParameter('subtractBaseline',true,@islogical);
p.addParameter('iampColor',[0,0,0],@isvector);
p.addParameter('indivBootCRF',[],@ismatrix);
p.parse(analysisParams,crfPlotStruct,crfStimulus,iampPoints,iampError,varargin{:});

if ~isempty(p.Results.indivBootCRF)
    indivBootCRF = p.Results.indivBootCRF;
end

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
            iampVals = iampPoints.paramMainMatrix(1:iampIndx)';
            if isfield(theModelResp, 'shaddedErrorBars')
                shdErrVals = theModelResp.shaddedErrorBars(:,1:modelIndx);
            end
            if exist('iampError','var')
                errVals = iampError.paramMainMatrix(:,1:iampIndx)';
            end
            if ~isempty(p.Results.indivBootCRF)
                theBoots = indivBootCRF(:,1:modelIndx)';
            end
        else
            crfValues = theModelResp.values((ii-1)*modelIndx+1:ii*modelIndx);
            xAxisModels = contrastSpacing((ii-1)*modelIndx+1:ii*modelIndx);
            iampVals = iampPoints.paramMainMatrix((ii-1)*iampIndx+1:ii*iampIndx)';
            if isfield(theModelResp, 'shaddedErrorBars')
                shdErrVals = theModelResp.shaddedErrorBars(:,(ii-1)*modelIndx+1:ii*modelIndx);
            end
            if exist('iampError','var')
                errVals = iampError.paramMainMatrix(:,(ii-1)*iampIndx+1:ii*iampIndx)';
            end
            if ~isempty(p.Results.indivBootCRF)
                theBoots = indivBootCRF(:,(ii-1)*modelIndx+1:ii*modelIndx)';
            end
        end
        if p.Results.subtractBaseline
            offestVal = iampPoints.paramMainMatrix(end);
            crfValues = crfValues - offestVal;
            iampVals  = iampVals - offestVal;
            if ~isempty(p.Results.indivBootCRF)
                theBoots  = theBoots - offestVal;
            end
        end
        xAxisIamp = maxConVal.*analysisParams.contrastCoding;
        %% Plot the stuff
        subplot(rws,cols,ii); hold on
        if ~isempty(p.Results.indivBootCRF)
            plot(xAxisModels,theBoots,'--','Color',[0.3, 0.3, 0.3,.4],'LineWidth', 0.75);
        end
        h(jj) = plot(xAxisModels,crfValues,'color',theModelResp.plotColor,'LineWidth', 1.0);
        if isfield(theModelResp, 'shaddedErrorBars')
            shadedErrorBars(xAxisModels,crfValues,shdErrVals,'lineprops',{'color',theModelResp.plotColor});
        end
        
        if exist('iampError','var')
            if size(errVals,2) == 2
                upperCI = errVals(:,1);
                lowerCI = errVals(:,2);
                q1    = errorbar(xAxisIamp, iampVals, lowerCI, upperCI, 'o','MarkerSize',7,...
                    'MarkerEdgeColor','k', 'MarkerFaceColor', p.Results.iampColor, ...
                    'LineWidth',1.0,'Color','k');
            else
                                q1    = errorbar(xAxisIamp,iampVals, errVals, 'o','MarkerSize',7,...
                    'MarkerEdgeColor','k', 'MarkerFaceColor', p.Results.iampColor, ...
                    'LineWidth',1.0,'Color','k');
            end
        else
            q1    = scatter(xAxisIamp,iampVals, 36,'o','MarkerFaceColor', p.Results.iampColor, ...
                'MarkerEdgeColor', 'k');
        end
        % put info
        ylabel('Mean Beta Weight')
        xlabel('Contrast')
        title(sprintf('LM stim = %s', num2str(analysisParams.LMVectorAngles(ii))));
        ylim([-0.3 1.4]);
        xlim([0 0.6]);
        set(gca,'xscale','log')
        set(gca, 'FontName', 'Helvetica', 'FontSize', 14,'FontWeight', 'normal');
        
    end
end

fields{end+1}  = 'IAMP Points';
legend([h, q1], fields, 'Location','NorthWest')
set(gcf, 'Position',  [0, 0, 1800, 1300])
end


