function plotTimeCourse(analysisParams, )

% Provides a higher resolution contrast/direcetions base for CRF predictions
% 
% Syntax:
%   [crfStimulus] = upsampleCRF(analysisParams)
%             
% Description:
%   This function takes in a fitting object, parameters, and a cell of packets
%   and returns the timecourse prediction of the model. 
%
% Inputs:    
%   analysisParams      - Struct of important information for the
%                         analysis. Relevant fields for this are:
%                           * contrastCoding 
%                           * directionCoding 
%                           * maxContrastPerDirection 
%                           * theDimention 
%                           * numSamples - upsample resolution 
% Outputs:
%   crfStimulus         - Upsampled contrast/directions stimuli. 
%
% Optional key/value pairs:
%   none

% History:
%   01/23/2019 MAB Wrote it. 

% subplot size
rws = ceil(sqrt(length(packets)));
cols = rws-1;
if rws*cols < length(packets)
    cols = rws;
end

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


legend([p, q1], fields, 'IAMP Points')
end


