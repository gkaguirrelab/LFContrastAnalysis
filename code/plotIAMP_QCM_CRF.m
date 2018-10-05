function [] = plotIAMP_QCM_CRF(analysisParams,meanIAMPBetas,semIAMPBetas,paramsQCMFit)
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

% Generate prediction to stimuli based on QCM fit to arbitrary stim

contrastSpacing = linspace(max(analysisParams.contrastCoding),min(analysisParams.contrastCoding),analysisParams.numSamples);
QCMStim.values = [generateStimCombinations(contrastSpacing,analysisParams.directionCoding,analysisParams.maxContrastPerDir,analysisParams.theDimension)];
QCMStim.timebase = linspace(1,max(length(meanIAMPBetas(1:end-1))),length(QCMStim.values));
temporalFitQCM = tfeQCM('verbosity','none','dimension',analysisParams.theDimension);
QCMResponses = computeResponse(temporalFitQCM,paramsQCMFit,QCMStim,[]);

% Subplot size
rws = ceil(sqrt(size(analysisParams.directionCoding,2)));
cols = rws;
indx = length(analysisParams.contrastCoding);
figure
for ii = 1:size(analysisParams.directionCoding,2)
    
    aa = (length(QCMResponses.values))/size(analysisParams.directionCoding,2);
    if ii == 1
        betas = meanIAMPBetas(1:indx);
        error = semIAMPBetas(1:indx);
        qcmSmooth = QCMResponses.values(1:aa);
    else
        betas = meanIAMPBetas((ii-1)*indx+1:ii*indx);
        error = semIAMPBetas((ii-1)*indx+1:ii*indx);
        qcmSmooth = QCMResponses.values((ii-1)*aa+1:ii*aa);
    end
    
    subplot(rws,cols,ii); hold on
    errorbar(analysisParams.contrastCoding*100,betas,error)
    plot(contrastSpacing*100,qcmSmooth)
    
    if isequal(analysisParams.directionCoding(:,ii),[1;1;0])
        pos = find(ismember(analysisParams.directionCoding',[1,1,0],'rows'));
        title(sprintf('L+M: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    elseif isequal(analysisParams.directionCoding(:,ii),[1;-1;0])
        pos = find(ismember(analysisParams.directionCoding',[1,-1,0],'rows'));
        title(sprintf('L-M: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    elseif isequal(analysisParams.directionCoding(:,ii),[1;0;0])
        pos = find(ismember(analysisParams.directionCoding',[1,0,0],'rows'));
        title(sprintf('L Isolating: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    elseif isequal(analysisParams.directionCoding(:,ii),[0;1;0])
        pos = find(ismember(analysisParams.directionCoding',[0,1,0],'rows'));
        title(sprintf('M Isolating: Max Contrast = %s',num2str(analysisParams.maxContrastPerDir(pos))))
    else
        title(sprintf('%s Degrees: Max Contrast = %s',num2str(analysisParams.LMVectorAngles(ii)),num2str(analysisParams.maxContrastPerDir(ii))))
    end
    ylabel('Mean Beta Weight')
    xlabel('Percent of Max Contrast')
    ylim([-0.2 1]);
end

end