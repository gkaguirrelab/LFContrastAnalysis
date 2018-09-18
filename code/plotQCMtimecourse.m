function [] = plotQCMtimecourse(paramsFitIAMP,packetPocket,meanIAMPBetas,analysisParams,fitResponseStructQCM);
% Takes in a text file name and retuns a cell of the lines of the text file
%
% Syntax:
%   filesCell = textFile2cell(inFile)
%
% Description:
%    This function takes in a file name for a trext file and returns a cell
%    that is composed of the lines of the text file. Example of this would
%    be a text file of file names the output is a cell of files names.
%
% Inputs:
%    inFile            - File name of a text file. (string)
%
% Outputs:
%    fileCell          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    none

% MAB 09/09/18

% IAMP object
temporalFitIAMP = tfeIAMP('verbosity','none');

% Get subplot sizing 
rws = ceil(sqrt(analysisParams.numAcquisitions));
cols = rws-1;
if rws*cols < analysisParams.numAcquisitions
    cols = rws;
end

% Open figure
figure

% Use QCM fit to IAMP to predict timecourse.
for jj = 1:analysisParams.numAcquisitions   
    % Regenerate IAMP predictions to time series
    IAMPResponses = temporalFitIAMP.computeResponse(paramsFitIAMP{jj},packetPocket{jj}.stimulus,packetPocket{jj}.kernel);
    
    % Plot them
    subplot(rws,cols,jj); hold on
    plot(packetPocket{jj}.response.timebase,packetPocket{jj}.response.values,'Color',[1 0 0]);
    plot(IAMPResponses.timebase, IAMPResponses.values,'Color',[0 1 0]);
    
    % Doctor up the parameters to use mean IAMP values and plot again
    paramsFitIAMPMean = paramsFitIAMP{jj};
    paramsFitIAMPMean.paramMainMatrix(1:21) = meanIAMPBetas(1:end-1); %end minus 1 for the attentional events
    IAMPResponsesMean = temporalFitIAMP.computeResponse(paramsFitIAMPMean,packetPocket{jj}.stimulus,packetPocket{jj}.kernel);
    plot(IAMPResponsesMean.timebase,IAMPResponsesMean.values,'Color',[0 0.5 1]);
   
    % Doctor up parameters to use the QCM fit to the mean IAMP
    paramsFitIAMPQCM = paramsFitIAMP{jj};
    paramsFitIAMPQCM.paramMainMatrix(1:21) = fitResponseStructQCM.values';
    IAMPResponsesQCM = temporalFitIAMP.computeResponse(paramsFitIAMPQCM,packetPocket{jj}.stimulus,packetPocket{jj}.kernel);
    plot(IAMPResponsesQCM.timebase,IAMPResponsesQCM.values,'Color',[0 0 0]);
    
    % Set axis labels
    ylabel('PSC')
    xlabel('Time (mS)')
    
    % Change line size
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
end

end