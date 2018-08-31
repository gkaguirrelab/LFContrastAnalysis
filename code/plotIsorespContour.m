function [] = plotIsorespContour(QCMparams,IAMPBetas,contrastLevels,directionCoding,thresh)
% Plots an isorepsonse contour for a given 2D ellipse fit along the data points
%
% Syntax:
%   [] = plotIsorespContour(varargin)
%
% Description:
%    Plots an isoresponse contour using the QCM fits and IAMP data points
%
% Inputs:
%    QCMParams       - paramsFit outpus from the QCM fit response function
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

%% Inerpolate the IAMP CRF
for ii = 1:length(IAMPBetas)
    contrast = interp1(IAMPBetas{ii},contrastLevels{ii},thresh,'pchip');
    dataPoints(ii,1:2) = contrast.*directionCoding{ii};
end

%% Plot data points 
figure; hold on
sz = 50;
scatter(dataPoints(:,1),dataPoints(:,2),sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ylim([-1, 1])
xlim([-1, 1])
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = ':'; % dotted
line(get(axh,'XLim'), [0 0], 'Color', color, 'LineStyle', linestyle);
line([0 0], get(axh,'YLim'), 'Color', color, 'LineStyle', linestyle);
xlabel('L Contrast')
ylabel('M Contrast')

end