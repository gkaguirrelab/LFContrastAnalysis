function figHdl = plotTimeCourse(analysisParams,timeCoursePlot, baselineShift, numSubPlots)
%Plots the predictions of each model for each run.
%
% Syntax:
%   plotTimeCourse(analysisParams,timeCoursePlot, baselineShift)
%
% Description:
%   Plots the predictions of each model for each run.
%
% Inputs:
%   analysisParams      - Struct of important information for the
%                         analysis.
%   timeCoursePlot      - A stuct of model predictions. Each field should
%                         be a different model containing subfields for
%                         avalues, timebase, and plotColor.
%   baselineShift       - A matrix with the shift needed to add the baseline
%                         back to the time course predicitions. This should be a
%                         matrix of size numSessions x numAcquisitions. If
%                         no baseline shisft have a matrix of all zeros
%   numSubPlots         - Number of subplots
% Outputs:
%   figHdl              - Figure handle
%
% Optional key/value pairs:
%   none

% History:
%   01/24/2019 MAB Wrote it.

% subplot size
rws = ceil(sqrt(numSubPlots));
cols = rws-1;
if rws*cols < numSubPlots
    cols = rws;
end

% reshape baselineShift
baselineShift = baselineShift';
baselineShift = baselineShift(:);

fields = fieldnames(timeCoursePlot);

figHdl = figure;

for ii = 1:numSubPlots
    
    for jj = 1:length(fields)
        
        % get the sturct
        theModelResp = eval(['timeCoursePlot.', fields{jj}]);
        
        % get the time course
        if strcmp('timecourse',fields{jj})
            response = theModelResp{ii}.values;
            timeCourse = response;
        elseif strcmp('qcm',fields{jj})
            response = theModelResp{ii}.values + baselineShift(ii);
            qcmPred = response;
        elseif strcmp('IAMP',fields{jj})
            response = theModelResp{ii}.values + baselineShift(ii);
            iampResp = response;
        else
            response = theModelResp{ii}.values + baselineShift(ii);
        end
        
        subplot(rws,cols,ii); hold on
        if isfield(theModelResp{ii}, 'shaddedErrorBars')
            shadedErrorBars(theModelResp{ii}.timebase,response,theModelResp{ii}.shaddedErrorBars,...
                'lineprops',{'color',theModelResp{ii}.plotColor},'patchSaturation', 0.25);
        end
        p(jj) = plot(theModelResp{ii}.timebase,response,'color',theModelResp{ii}.plotColor,'LineWidth', 1.0);
        
        
        
    end
    
    if ~isempty(find(strcmp(fields(:), 'qcm')))
        tempR2 = corr([timeCourse',qcmPred']).^2;
        qcmR2 = tempR2(1,2);
    end
    if ~isempty(find(strcmp(fields(:), 'IAMP')))
        tempR2 = corr([timeCourse',iampResp']).^2;
        iampR2 = tempR2(1,2);
    end
    
    % put info
    ylabel('PSC')
    xlabel('Time mS')
    if ~isempty(find(strcmp(fields(:), 'IAMP')))
        title(sprintf('Run %s: QCM R^{2} %s, IAMP R^{2} %s', num2str(ii),num2str(round(qcmR2,2)),num2str(round(iampR2,2))));
    else
        title(sprintf('Run %s', num2str(ii)));
    end
    set(gca, 'FontName', 'Helvetica', 'FontSize', 14,'FontWeight', 'normal');
end


legend(p, fields)
set(gcf, 'Position',  [0, 0, 1800, 1300])
end


