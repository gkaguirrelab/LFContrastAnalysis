function [hdl] = plotIsoresponse(analysisParams,meanIAMPBetas,paramsQCMFit,thresholds,nrParams,colors)
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
        sortedBetas{ii} = meanIAMPBetas.paramMainMatrix(1:numContrast);
    else
        sortedBetas{ii} = meanIAMPBetas.paramMainMatrix((ii-1)*numContrast+1:ii*numContrast);
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

% Dummy vars for plotting
hdl = [];
hleglines = [];

% Loop over thesholds and plot the ellipse fits on the same plot
for jj = 1:length(thresholds)
    color = colors(jj,:);
    [hdl,scatterHdl] = plotIsorespContour(paramsQCMFit,nrParams,sortedBetas,analysisParams, directions,thresholds(jj),hdl,color);
    hleglines = [hleglines scatterHdl];
    legendNames{jj} = num2str(thresholds(jj));
end

xlim([-.5 .5])
ylim([-.5 .5])
axh = gca; % use current axes
axisColor = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', axisColor, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', axisColor, 'LineStyle', linestyle);
xlabel('L Contrast')
ylabel('M Contrast')
legend(hleglines,legendNames)
title('Isoresponse Contour')
axis square
set(gcf, 'Position',  [0, 0, 900, 900])


end