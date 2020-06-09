subjIds = {'KAS25','KAS25_replication','AP26','AP26_replication','LZ23','LZ23_replication'};
analysisParams = getSubjectParams('KAS25');

for ii = 1:length(subjIds)
    [regLineParams] =  analyzeLFContrast_GLM_HighLow(subjIds{ii});
    
    slopes_lMinusM_high(ii) = regLineParams.high.lMinusM(2,1);
    offsets_lMinusM_high(ii)= regLineParams.high.lMinusM(1,1);
    
    slopes_lPlusM_high(ii) = regLineParams.high.lPlusM(2,1);
    offsets_lPlusM_high(ii)= regLineParams.high.lMinusM(1,1);
    
    slopes_lMinusM_low(ii) = regLineParams.low.lMinusM(2,1);
    offsets_lMinusM_low(ii)= regLineParams.low.lMinusM(1,1);
    
    slopes_lPlusM_low(ii) = regLineParams.low.lPlusM(2,1);
    offsets_lPlusM_low(ii)= regLineParams.low.lPlusM(1,1);
    
end

meanSlope_lMinusM_high  = mean(slopes_lMinusM_high);
meanOffset_lMinusM_high = mean(offsets_lMinusM_high);

meanSlope_lPlusM_high   = mean(slopes_lPlusM_high);
meanOffset_lPlusM_high  = mean(offsets_lPlusM_high);

meanSlope_lMinusM_low   = mean(slopes_lMinusM_low);
meanOffset_lMinusM_low  = mean(offsets_lMinusM_low);

meanSlope_lPlusM_low    = mean(slopes_lPlusM_low);
meanOffset_lPlusM_low   = mean(offsets_lPlusM_low);

regLine_lMinusM_high = @(x) meanSlope_lMinusM_high.*x + meanOffset_lMinusM_high;
regLine_lPlusM_high  = @(x) meanSlope_lPlusM_high.*x + meanOffset_lPlusM_high;
regLine_lMinusM_low  = @(x) meanSlope_lMinusM_low.*x + meanOffset_lMinusM_low;
regLine_lPlusM_low   = @(x) meanSlope_lPlusM_low.*x + meanOffset_lPlusM_low;

xPts = 0:1:20;

yPts_lMinusM_high = regLine_lMinusM_high(xPts);
yPts_lPlusM_high  = regLine_lPlusM_high(xPts);
yPts_lMinusM_low  = regLine_lMinusM_low(xPts);
yPts_lPlusM_low   = regLine_lPlusM_low(xPts);

%% PLOT IT

% plotting color

lMinusMColorLine = [244,194,194]./255;
lPlusMColorLine  = [204,204,255]./255;

figHndl = figure;

% Set the figure's size in inches
figureSizeInches = [20 12];
figHndl.Units = 'inches';


subplot(1,2,1);
hold on;
axis square

l1 = line(xPts,yPts_lMinusM_high,'Color',lMinusMColorLine,'LineWidth',3);
l2 = line(xPts,yPts_lPlusM_high ,'Color',lPlusMColorLine,'LineWidth',3);

ylim([-0.5,1.0])
xlabel('Eccentricity')
ylabel('GLM Beta Weight')
title('High Contrast Condition');
set(gca, ...
    'XColor', [0.2 0.2 0.2], ...
    'YColor', [0.2 0.2 0.2], ...
    'FontName', 'Helvetica', ...
    'FontSize', 14, ...
    'FontWeight', 'normal', ...
    'TickLength',[0.01 0.01], ...
    'TickDir', 'out', ...
    'LineWidth', 0.7, ...
    'Box', 'off');


subplot(1,2,2);
axis square
hold on;
l3 = line(xPts,yPts_lMinusM_low,'Color',lMinusMColorLine,'LineWidth',3);
l4 = line(xPts,yPts_lPlusM_low ,'Color',lPlusMColorLine,'LineWidth',3);

ylim([-0.5,1.0])
xlabel('Eccentricity')
ylabel('GLM Beta Weight')
title('Low Contrast Condition');
set(gca, ...
    'XColor', [0.2 0.2 0.2], ...
    'YColor', [0.2 0.2 0.2], ...
    'FontName', 'Helvetica', ...
    'FontSize', 14, ...
    'FontWeight', 'normal', ...
    'TickLength',[0.01 0.01], ...
    'TickDir', 'out', ...
    'LineWidth', 0.7, ...
    'Box', 'off');
legend([l3,l4],{'L-M', 'L+M'})

set(figHndl, 'PaperSize',figureSizeInches);
set(figHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
% Full file name
figName =  fullfile(getpref(analysisParams.projectName,'figureSavePath'), ...
    [analysisParams.expSubjID,'_Average_scatter_High_Low_Contrast_GLM_hcp.pdf']);
% Save it
print(figHndl, figName, '-dpdf', '-r300');