function [hdl,scatterHdl] = plotIsorespContour(paramsQCM,nrParams,IAMPBetas,analysisParams, directionCoding,thresh,hdl,color)
% Plots an isorepsonse contour for a given 2D ellipse fit along the data points
%
% Syntax:
%   [] = plotIsorespContour(varargin)
%
% Description:
%    Plots an isoresponse contour using the QCM fits and IAMP data points
%
% Inputs:
%    paramsQCM       - paramsFit outpus from the QCM fit response function
%    IAMPBetas       - beta values for each
%    contrastLevels  - contrast values corresponding to each beta weight in each direction
%    directionCoding - coding for directions in the XY plane e.g. [1,1] = L+M
%
% Outputs:
%    myRes          - The calculated difference between the two
%                     provided integer values.
%
% Optional key/value pairs:
%    None.
%

%% Parameters
% number of points for the ellipse
nQCMPoints = 100;

% Chose random colors if not defined
if isempty(color)
    color = [1,1,1];
    while sum(color) > 2.4
        color = rand(1,3);
    end
end

%% Inerpolate the IAMP CRF using the  naka rushton fits to find the contrast value that corresponds with the threshold
for ii = 1:size(nrParams,1)
    
    % Invert Naka-Rushton function to get the contrast value that
    %  Rmax  = params(1)
    %  sigma = params(2)
    %  n     = params(3)
    maxConVal = analysisParams.maxContrastPerDir(ii);
    maxContrastSpacing = maxConVal.*analysisParams.contrastCoding;
    if thresh <= nrParams(1)
        contrastsNR(ii) = InvertNakaRushton([nrParams(ii,1),nrParams(ii,2),nrParams(ii,3)],thresh);
        contrastsLI(ii) = interp1(IAMPBetas{ii},maxContrastSpacing',thresh);
    else
        contrastsNR(ii) = NaN;
        contrastsLI(ii) = interp1(IAMPBetas{ii},maxContrastSpacing',thresh);
    end
      
    % Get the L,M plane coordinates by mulitplying the contrast needed by the direction coding.
    % NOTE: MB: I think this should be the sin and cos comp. of the
    % direction and not the coding.
    dataPointsNR(ii,1:2) = contrastsNR(ii).*directionCoding{ii};
    dataPointsLI(ii,1:2) = contrastsLI(ii).*directionCoding{ii};
end

%% Compute QCM ellipse to the plot
%
% Step 1. Invert Naka-Rushton to go from thresh back to
% corresponding equivalent contrast.
eqContrast = InvertNakaRushton([paramsQCM.crfAmp,paramsQCM.crfSemi,paramsQCM.crfExponent],thresh);
circlePoints = eqContrast*UnitCircleGenerate(nQCMPoints);
[~,Ainv,Q] = EllipsoidMatricesGenerate([1 paramsQCM.Qvec],'dimension',2);
ellipsePoints = Ainv*circlePoints;
checkThresh = ComputeNakaRushton([paramsQCM.crfAmp,paramsQCM.crfSemi,paramsQCM.crfExponent],diag(sqrt(ellipsePoints'*Q*ellipsePoints)));
if (any(abs(checkThresh-thresh) > 1e-10))
    error('Did not invert QCM model correctly');
end

%% Plot data points
if (isempty(hdl))
    hdl = figure; hold on
else
    figure(hdl); hold on
end
sz = 50;
scatterHdl = scatter(dataPointsNR(:,1),dataPointsNR(:,2),sz,'MarkerEdgeColor',color,'MarkerFaceColor',color,'LineWidth',1.5);
scatter(dataPointsLI(:,1),dataPointsLI(:,2),sz,color,'x')


% Add ellipse
plot(ellipsePoints(1,:),ellipsePoints(2,:),'color', color);


end
