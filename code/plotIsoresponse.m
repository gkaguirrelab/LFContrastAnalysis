function [hdl] = plotIsoresponse(analysisParams,meanIAMPBetas,paramsQCMFit,thresholds)
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

% Get the number of conditions (directions) 
numCond = size(analysisParams.directionCoding,2);

% Get the number of contrast levels
numContrast = length(analysisParams.contrastCoding );


% Extract the beta weights for each direction (one weight per contrast
% level and condtion)
for ii = 1:numCond
    if ii == 1
        sortedBetas{ii} = meanIAMPBetas(1:numContrast);
    else
        sortedBetas{ii} = meanIAMPBetas((ii-1)*numContrast+1:ii*numContrast);
    end
    
    % Scale the contrast spacing by the maximum contrast per direction 
    contrasts{ii} = analysisParams.contrastCoding*analysisParams.maxContrastPerDir(ii);
end

% Get the direction coding and take only the L and M coding  
directions = analysisParams.directionCoding;
if analysisParams.theDimension == 2 & size(analysisParams.directionCoding,1) > 2
    directions(3,:) = [];
end

% Turn the direction coding matrix in cell array 
directions = mat2cell(directions,size(directions,1),ones(1,size(directions,2)));

hdl = [];
hleglines = []
for jj = 1:length(thresholds)
    [hdl,scatterHdl] = plotIsorespContour(paramsQCMFit,sortedBetas,contrasts,directions,thresholds(jj),hdl,[]);
    hleglines = [hleglines scatterHdl];
    legendNames{jj} = num2str(thresholds(jj));
end
legend(hleglines,legendNames)
title('Isoresponse Contour')
end