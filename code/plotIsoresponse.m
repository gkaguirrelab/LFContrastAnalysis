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

numCond = size(analysisParams.directionCoding,2);
numContrast = length(analysisParams.contrastCoding );

for ii = 1:numCond
    if ii == 1
        sortedBetas{ii} = meanIAMPBetas(1:numContrast) - meanIAMPBetas(end-1);
    else
        sortedBetas{ii} = meanIAMPBetas((ii-1)*numContrast+1:ii*numContrast) - meanIAMPBetas(end-1);
    end
    
    contrasts{ii} = analysisParams.contrastCoding*analysisParams.maxContrastPerDir(ii)
end

directions = analysisParams.directionCoding;
if analysisParams.theDimension == 2 & size(analysisParams.directionCoding,1) > 2
    directions(3,:) = [];
end

directions = mat2cell(directions,size(directions,1),ones(1,size(directions,2)));

hdl = [];
for jj = 1:length(thresholds)
    hdl = plotIsorespContour(paramsQCMFit,sortedBetas,contrasts,directions,thresholds(jj),hdl,[]);
end

end