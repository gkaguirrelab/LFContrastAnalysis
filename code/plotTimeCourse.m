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
        
        theModelResp = eval(['timeCoursePlot.', fields{jj}]);
        
        
        if strcmp('timecourse',fields{jj})
            response = theModelResp{ii}.values;
        else
            response = theModelResp{ii}.values + baselineShift(ii);
            
        end
        
        subplot(rws,cols,ii); hold on
        p(jj) = plot(theModelResp{ii}.timebase,response,'color',theModelResp{ii}.plotColor);
        
        
        % put info
        ylabel('PSC')
        xlabel('Time mS')
        title(sprintf('Run = %s', num2str(ii)));
        
        
    end
end


legend(p, fields)
set(gcf, 'Position',  [0, 0, 1800, 1300])
end


