function [figHdl] = plotCRF(analysisParams, crfPlotStruct, crfStimulus, iampPoints, iampSEM, varargin)
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
%    iampSEM                   - Error bars for iampPoints
%
% Outputs:
%    figHdl                    - Figure handle
%
% Optional key/value pairs:
%    subtractBaseline          - The baseline value to be subtracted from the CRF values

% MAB 09/09/18

% Subplot size
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addRequired('crfPlotStruct',@isstruct);
p.addRequired('crfStimulus',@isstruct);
p.addRequired('iampPoints',@isstruct);
p.addRequired('iampSEM',@isstruct);
p.addParameter('subtractBaseline',true,@islogical);


p.parse(analysisParams,crfPlotStruct,crfStimulus,iampPoints,iampSEM,varargin{:});



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
            if exist('iampSEM','var')
                errVals = iampSEM.paramMainMatrix(1:iampIndx)';
            end
        else
            crfValues = theModelResp.values((ii-1)*modelIndx+1:ii*modelIndx);
            xAxisModels = contrastSpacing((ii-1)*modelIndx+1:ii*modelIndx);
            iampVals = iampPoints.paramMainMatrix((ii-1)*iampIndx+1:ii*iampIndx)';
            if exist('iampSEM','var')
                errVals = iampSEM.paramMainMatrix((ii-1)*iampIndx+1:ii*iampIndx)';
            end
        end
        if p.Results.subtractBaseline
            %offestVal = crfValues(end);
            offestVal = iampPoints.paramMainMatrix(end);
            crfValues = crfValues - offestVal;
            iampVals = iampVals - offestVal;
        end
        xAxisIamp = maxConVal.*analysisParams.contrastCoding;
        %% Plot the stuff
        subplot(rws,cols,ii); hold on
        h(jj) = plot(xAxisModels,crfValues,'color',theModelResp.color);
        q1    = scatter(xAxisIamp,iampVals,'ok');
        if exist('iampSEM','var')
            errorbar(xAxisIamp,iampVals, errVals, 'LineStyle','none');
        end
        % put info
        ylabel('Mean Beta Weight')
        xlabel('Contrast')
        title(sprintf('LM stim = %s', num2str(analysisParams.LMVectorAngles(ii))));
        ylim([-0.3 1.4]);
        
    end
end

fields{end+1}  = 'IAMP Points';
legend([h, q1], fields, 'Location','NorthWest')
set(gcf, 'Position',  [0, 0, 1800, 1300])
end


