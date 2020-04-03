function t_MichaelFigure
    % Get the data to plot
    [theEllipse, xSampleBase, nrVals, ...
        eqContrastPts, semQCMParams, qcmParams] = getData();
    
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();
    
    % Assemble the color order for all sequential data sets
    theColorOrder = eqContrastPts.colorVals;
    
    % Line color
    theLineColor = [0.4078, 0.2784, 0.5765];
    
    % Apply the default recipe with some overrides
    plotlabOBJ.applyRecipe(...
        'colorOrder', theColorOrder, ...
        'lineMarkerSize', 7.5, ...
        'scatterMarkerFaceAlpha', 0.6, ...
        'axesTickDir', 'in', ...
        'axesTickLength', [.02 .02], ...
        'axesLineWidth', 2.0, ...
        'axesTitleFontWeight', 'bold', ...
        'axesFontSize', 15, ...
        'axesLabelFontSizeMultiplier', 1.1, ...
        'figureWidthInches',11, ...
        'figureHeightInches', 5);
    
    % New figure
    hFig = figure(1); clf;
    
    % Generate axes in a [1x2] layout
    theAxesGrid = plotlab.axesGrid(hFig, ...
        'rowsNum', 1, 'colsNum', 2, ...
        'rightMargin', 0.15, ...
        'leftMargin', 0.08, ...
        'widthMargin', 0.08, ...
        'bottomMargin', 0.08, ...
        'topMargin', 0.05);
    
    % -------- The left plot (1,1) --------
    % Retrieve the axes
    theCurrentAxes = theAxesGrid{1,1};
    
    % Plot the cross hairs
    plotlab.crossHairs2D(theCurrentAxes, ...
        'xRange', [-1 1], 'yRange', [-1 1], ...
        'LineWidth', 2, 'LineStyle', ':', ...
        'LineColor', [.3 .3 .3]);
    
    % Plot the ellipse
    line(theCurrentAxes,theEllipse.x,theEllipse.y, 'Color', theLineColor)
    
    % Text containing math set in LaTeX
    if isempty(semQCMParams)
        modelTxtTheta = ['${\theta} = ' num2str(round(qcmParams.Qvec(2),2)) '^{\circ}$'];
        modelTxtMAR   = ['$m_ratio = ' num2str(round(qcmParams.Qvec(1),2)) '$'];
    else
        modelTxtTheta = ['${\theta} = ' num2str(round(qcmParams.Qvec(2),2)) '^{\circ} {\pm} ' ...
            num2str(round(semQCMParams.Qvec(2),2)) '$'];
        modelTxtMAR   = ['$m = ' num2str(round(qcmParams.Qvec(1),2)) '{\pm} ' ...
            num2str(round(semQCMParams.Qvec(1),2)) '$'];
    end
    
    % Add the above text to the plot
    theTextHandle = text(theCurrentAxes, -.9,.9 , modelTxtTheta, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);
    theTextHandle = text(theCurrentAxes, -.9,.76 , modelTxtMAR, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);
    
    % Axes adjustments
    axis(theCurrentAxes,'square');
    set(theCurrentAxes, 'XLim', [-1 1], 'YLim', [-1 1], ...
        'XTick', -1:.5:1, 'YTick', -1:.2:1, 'YGrid', 'off');

    
    % Labels
    xlabel(theCurrentAxes, 'L contrast');
    ylabel(theCurrentAxes,'M contrast');
    
    % Title
    title(theCurrentAxes, 'Isoresponse Contour');
    % --------------------------------------
    
    
    % -------- The right plot (1,2) --------
    % Retrieve the axes
    theCurrentAxes = theAxesGrid{1,2};
    
    hold(theCurrentAxes, 'on');
    % plot each point with its associated color
    for ii = 1:length(eqContrastPts.values)
        % plot it
        scatter(theCurrentAxes, eqContrastPts.values(1,ii),eqContrastPts.values(2,ii),...
            'LineWidth', 1.0);
    end
    
    % The non-linearity
    plot(theCurrentAxes, xSampleBase, nrVals, '-', 'Color', theLineColor);
    
    % Set color map range
    caxis(theCurrentAxes, [-45,112.5]);
    
    % Set the color map to custom 8 color
    colormap(theCurrentAxes, eqContrastPts.colorMap);
    
    % Add colorbar without resizing the figure
    plotlab.colorbar(theCurrentAxes,'Chromatic Direction (angles in L/M plane)', 12, ...
        'Location','eastoutside' ,...
        'Ticks',[-35,-15,5,25,42.5,62.5,82.5,102.5],...
        'TickLabels',{'-45^o','-22.5^o','0^o','22.5^o','45^o','67.5^o','90^o','112.5^o'});
    
    % Text containing math set in LaTeX
    if isempty(semQCMParams)
        modelTxtAmp  = ['$Amp = ' num2str(round(qcmParams.crfAmp,2)) '$'];
        modelTxtExp  = ['$Exp = ' num2str(round(qcmParams.crfExponent,2)) '$'];
        modelTxtSemi = ['$Semi = ' num2str(round(qcmParams.crfSemi,2)) '$'];
    else
        modelTxtAmp  = ['$Amp = ' num2str(round(qcmParams.crfAmp,2)) '{\pm} ' ...
            num2str(round(semQCMParams.crfAmp,2)) '$'];
        modelTxtExp  = ['$Exp = ' num2str(round(qcmParams.crfExponent,2)) '{\pm} ' ...
            num2str(round(semQCMParams.crfExponent,2)) '$'];
        modelTxtSemi = ['$Semi = ' num2str(round(qcmParams.crfSemi,2)) '{\pm} ' ...
            num2str(round(semQCMParams.crfSemi,2)) '$'];
    end
    
    % Add the above text to the plot
    theTextHandle = text(theCurrentAxes, .0075,9.5 , modelTxtAmp, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);
    theTextHandle = text(theCurrentAxes, .0075,8.8, modelTxtExp, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);
    theTextHandle = text(theCurrentAxes, .0075,8.1 , modelTxtSemi, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 12, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);
    
    % Axes adjustments
    set(theCurrentAxes, 'XScale', 'log', ...
        'TickDir'     , 'out', ...
        'XLim', [0 1], 'YLim', [0 10], ...
        'XMinorTick'  , 'on', 'YMinorTick', 'on', ...
        'YTick'       , 0:2:10    , ...
        'XTick'       , [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0], ...
        'XTickLabel'  , {'0%','1%','2%','5%','10%','20%','50%','100%'});
    axis(theCurrentAxes,'square');
    
    % Labels
    xlabel(theCurrentAxes, 'Equivalent Contrast');
	ylabel(theCurrentAxes, 'Response');
    
    % Title
    title (theCurrentAxes, 'Response Nonlinearlity');
    
    % Export the figure to the gallery directory in PNG format
    plotlabOBJ.exportFig(hFig, 'png', 'twinPlots', 'gallery');
end

function [theEllipse, xSampleBase, nrVals, eqContrastPts, semQCMParams, qcmParams] = getData()
    % Load the data
    dataFileName = '/Users/michael/Desktop/vars4nicolas.mat';
    load(dataFileName,'ellipsePoints', 'eqContrastPts', ...
        'qcmParams', 'semQCMParams');
    
    % Subtract the baseline
    eqContrastPts.values(2,:) = eqContrastPts.values(2,:) - eqContrastPts.values(2,end);

    theEllipse.x = ellipsePoints(1,:);
    theEllipse.y = ellipsePoints(2,:);
    
    % Get the NR params
    amp  = qcmParams.crfAmp;
    exp  = qcmParams.crfExponent;
    semi = qcmParams.crfSemi;

    % Define NR function
    NR = @(c) amp .* ((c.^exp)./(c.^exp + semi.^exp));

    % set the sampling and range of the NR
    xSampleBase  = 0:0.01:1;

    % Evaluate the NR at the xSampleBase Points
    nrVals = NR(xSampleBase);
end
