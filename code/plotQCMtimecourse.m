function [] = plotQCMtimecourse(paramsFitIAMP,packetPocket,meanIAMPBetas,analysisParams,fitResponseStructQCM,baselineBetas);
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
rws = ceil(sqrt(analysisParams.numAcquisitions*length(analysisParams.sessionFolderName)));
cols = rws-1;
if rws*cols < analysisParams.numAcquisitions*length(analysisParams.sessionFolderName)
    cols = rws;
end

% Set indexing for betas 
betaLength = (length(meanIAMPBetas))/length(analysisParams.sessionFolderName);

% Open figure
figure

% Use QCM fit to IAMP to predict timecourse.
counter = 1;
for ii = 1:length(analysisParams.sessionFolderName)
    for jj = 1:analysisParams.numAcquisitions
        % Regenerate IAMP predictions to time series
        IAMPResponses = temporalFitIAMP.computeResponse(paramsFitIAMP{counter},packetPocket{counter}.stimulus,packetPocket{counter}.kernel);
        
        % Plot them
        subplot(rws,cols,counter); hold on
        plot(packetPocket{counter}.response.timebase,packetPocket{counter}.response.values,'Color',[1 0 0]);
        plot(IAMPResponses.timebase, IAMPResponses.values,'Color',[0 1 0]);
        
        % Doctor up the parameters to use mean IAMP values and plot again
        paramsFitIAMPMean = paramsFitIAMP{counter};
        paramsFitIAMPMean.paramMainMatrix(1:end-1) = [meanIAMPBetas(1+((ii-1)*betaLength):ii*betaLength);0] + baselineBetas(jj,ii); 
        IAMPResponsesMean = temporalFitIAMP.computeResponse(paramsFitIAMPMean,packetPocket{counter}.stimulus,packetPocket{counter}.kernel);
        plot(IAMPResponsesMean.timebase,IAMPResponsesMean.values,'Color',[0 0.5 1]);
        
        % Doctor up parameters to use the QCM fit to the mean IAMP
        paramsFitIAMPQCM = paramsFitIAMP{counter};
        paramsFitIAMPQCM.paramMainMatrix(1:end-1) = [fitResponseStructQCM.values(1+((ii-1)*betaLength):ii*betaLength), 0]' + baselineBetas(jj,ii);
        IAMPResponsesQCM = temporalFitIAMP.computeResponse(paramsFitIAMPQCM,packetPocket{counter}.stimulus,packetPocket{counter}.kernel);
        plot(IAMPResponsesQCM.timebase,IAMPResponsesQCM.values,'Color',[0 0 0]);
        
        % Set axis labels
        ylabel('PSC')
        xlabel('Time (mS)')
        title(sprintf('Session %s, Run %s', num2str(ii), num2str(jj)))
        % Change line size
        set(findall(gca, 'Type', 'Line'),'LineWidth',1);
        counter = counter +1; 
    end
end
legend('time course','IAMP fit',' Mean IAMP params', 'Mean QCM params')
end