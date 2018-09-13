function [] = plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas)
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

%% Generate prediction to stimuli based on QCM fit to arbitrary stim

contrastSpacing = linspace(max(analysisParams.contrastCoding),min(analysisParams.contrastCoding),analysisParams.numSamples);
QCMStim.values = [generateStimCombinations(contrastSpacing,directionCoding,maxContrastPerDir,theDimension),[0;0]];
QCMStim.timebase = linspace(1,max(thePacket.response.timebase),length(QCMStim.values));

QCMResponses = computeResponse(temporalFitQCM,paramsQCMFit,QCMStim,[]);

% Subplot size
rws = ceil(sqrt(size(analysisParams.directionCoding,2)));
cols = rws;
figure
for ii = 1:size(analysisParams.directionCoding,2)
    if ii == 1
        betas = meanIAMPBetas(1:indx)- meanIAMPBetas(end);
        error = semIAMPBetas(1:indx);
    else
        betas = meanIAMPBetas(indx+1:ii*indx) - meanIAMPBetas(end);
        error = semIAMPBetas(indx+1:ii*indx);
    end
    subplot(rws,cols,ii); hold on
    errorbar(xPos,betas,error)
    plot(contrastSpacing*100,QCMResponses.values(1:25)-QCMResponses.values(101))
    
    
    if isequal(analysisParams.directionCoding(:,ii),[1;1;0])
        [~,indx]=ismember(X,M,'rows')
        title(sprintf('L+M: Max Contrast = 6%')
    elseif isequal(analysisParams.directionCoding(:,ii),[1;-1;0])
        title('L-M: Max Contrast = 6%')
    elseif isequal(analysisParams.directionCoding(:,ii),[1;0;0])
        title('L-iso: Max Contrast = 6%')
    elseif isequal(analysisParams.directionCoding(:,ii),[0;1;0])
        title('M-iso: Max Contrast = 6%')
    else
    end
    ylabel('Mean Beta Weight')
    xlabel('Percent of Max Contrast')
    ylim([-0.2 1]);
end

end