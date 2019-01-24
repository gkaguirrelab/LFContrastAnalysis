function plotTimeCourse(analysisParams,timeCoursePlot, baselineShift)
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
% baselineShift         - A matrix with the shisft need to add the baseline
%                         back to the time course predicitions. should be a
%                         matrix of size numSessions x numAcquisitions. if 
%                         no baseline shisft have a matrix of all zeros 
% Outputs:
%   crfStimulus         - Upsampled contrast/directions stimuli. 
%
% Optional key/value pairs:
%   none

% History:
%   01/23/2019 MAB Wrote it. 

% subplot size
rws = ceil(sqrt(analysisParams.numSessions*analysisParams.numAcquisitions));
cols = rws-1;
if rws*cols < analysisParams.numSessions*analysisParams.numAcquisitions
    cols = rws;
end

% reshape baselineShift
baselineShift = baselineShift';
baselineShift = baselineShift(:);

fields = fieldnames(timeCoursePlot);

figHdl = figure; 

for ii = 1:analysisParams.numSessions*analysisParams.numAcquisitions
    
    for jj = 1:length(fields)
        
        theModelResp = eval(['timeCoursePlot.', fields{jj}]);
        
        response = theModelResp{ii}.values + baselineShift(ii);
        
        subplot(rws,cols,ii); hold on
        p(jj) = plot(theModelResp{ii}.timebase,response,'color',theModelResp{ii}.plotColor);

        
        % put info 
        ylabel('PSC')
        xlabel('Time mS')
        title(sprintf('Run = %s', num2str(ii)));


    end
end


legend(p, fields)
end


