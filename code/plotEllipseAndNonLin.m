function [figHndl] = plotEllipseAndNonLin(qcmParams, varargin)
% Plot the normalized ellipse and non-linearity from the QCM param
%
% Syntax:
%   [figHndl] = plotEllipseAndNonLin(qcmParams, varargin)
%
% Description:
%
%
% Inputs:
%    qcmParams         - Parameters from the QCM model fit.
%                           Angle, Minor Axis Ratio, NR Exp, NR AMP, NR
%                           Semi, NR Offset
%
% Outputs:
%    figHndl           - Figure Handle
%
% Optional key/value pairs:
%    nQCMPoints        - Resolution of the ellipse
%    ellipseColor      - 
% MAB 03/18/20

%% Input Parser
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('qcmParams',@isstruct);
p.addParameter('nQCMPoints',100,@isnumeric);
p.addParameter('plotColor',[0.3 0.45 0.81],@isvector);
p.addParameter('xSampleBase',[0:0.01:1],@isnumeric);
p.parse(qcmParams,varargin{:});

% Pull stuff out of the results struct
nQCMPoints   = p.Results.nQCMPoints;
plotColor = p.Results.plotColor;
xSampleBase  = p.Results.xSampleBase;

%% Ellipse Figure
% Get the eq. contrast needed for a normalized ellipse
desiredEqContrast = InvertNakaRushton([qcmParams.crfAmp,qcmParams.crfSemi,qcmParams.crfExponent],1);

% Generate a circle of calculated eq. contrast radius
circlePoints = desiredEqContrast*UnitCircleGenerate(nQCMPoints);

% Create transformation found from fitting the QCM
[~,Ainv,Q] = EllipsoidMatricesGenerate([1 qcmParams.Qvec],'dimension',2);

% Apply transformation
ellipsePoints = Ainv*circlePoints;

% Plot it
figHndl = figure;
subplot(1,2,1)
xlim([-1 1])
ylim([-1 1])

axh = gca; % use current axes
axisColor = [.3 .3 .3]; 
linestyle = ':'; % dotted
line([-1 1], [0 0], 'Color', axisColor, 'LineStyle', linestyle,'LineWidth', 2);
line([0 0], [-1 1], 'Color', axisColor, 'LineStyle', linestyle,'LineWidth', 2);
line(ellipsePoints(1,:),ellipsePoints(2,:),'color', plotColor, 'LineWidth', 2);
hXLabel = xlabel('L Contrast');
hYLabel = ylabel('M Contrast');
hTitle  = title('Isoresponse Contour');
set(gca,'FontSize',12);
set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 12);
set( hTitle, 'FontSize', 14,'FontWeight' , 'bold');
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -1:.2:1    , ...
    'LineWidth'   , 2         , ...
    'ActivePositionProperty', 'OuterPosition');
axis square


%% Non-linearity figure

% Get the NR params
amp  = qcmParams.crfAmp;
exp  = qcmParams.crfExponent;
semi = qcmParams.crfSemi;

% Define NR function
NR = @(c) amp .* ((c.^exp)./(c.^exp + semi.^exp));

% Evaluate the NR at the xSampleBase Points
nrVals = NR(xSampleBase);

% Plot it
subplot(1,2,2)
L1 = plot(xSampleBase,nrVals,'Color', plotColor, 'LineWidth', 2);

hTitle  = title ('Response Nonlinearlity');
hXLabel = xlabel('Equivalent Contrast'  );
hYLabel = ylabel('Response');

axis square

set(gca,'FontSize',12)
set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 12);
set( hTitle, 'FontSize', 14,'FontWeight' , 'bold');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:2:10    , ...
    'LineWidth'   , 2         , ...
    'ActivePositionProperty', 'OuterPosition');
ylim([0 10]);

set(gcf, 'Color', 'white' );


end