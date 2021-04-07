function plotHRF(subjId)
% plot the HRF for a subject
sessionNames = {subjId,[subjId, '_replication']};

for ii = 1:length(sessionNames)
    analysisParams = getSubjectParams(sessionNames{ii});
    analysisParams = loadHRF(analysisParams);
    
    theHrf{ii} = analysisParams.HRF;
    theHrf{ii}.timebase = theHrf{ii}.timebase/1000;
end

hrfFigHndl = figure; hold on

plot(theHrf{1}.timebase,theHrf{1}.values,'Color', [.4,.5,.1], 'LineWidth', 2);
plot(theHrf{2}.timebase,theHrf{2}.values,'Color', [.8,.55,.55], 'LineWidth', 2);

% set the axes and figure labels
hTitle  = title (sprintf('Hemodymanic Response Function - %s',subjId));
hXLabel = xlabel('Time (Seconds)');
hYLabel = ylabel('Amplitude');
set(gca,'FontSize',12)
set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 12);
set( hTitle, 'FontSize', 14,'FontWeight' , 'bold');

% format plot
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , [-1:.5:3] *1e-4    , ...
    'LineWidth'   , 2         , ...
    'ActivePositionProperty', 'OuterPosition');
ylim([-1 3] *1e-4);
xlim([0 32]);
set(gcf, 'Color', 'white' );
legend('Measurement Set 1', 'Measurement Set 2')

% save it
figSavePath = fullfile('/Users','mbarnett','labDropbox','LFContrastPaper','2019-LFContrastPaper','Figures','reviewFigs');

set(hrfFigHndl, 'Renderer', 'Painters');
figureSizeInches = [6.5 4];
set(hrfFigHndl, 'PaperUnits', 'inches');
set(hrfFigHndl, 'PaperSize',figureSizeInches);
set(hrfFigHndl, 'PaperPosition', [0 0 figureSizeInches(1) figureSizeInches(2)]);
figNameHrf = fullfile(figSavePath,[analysisParams.expSubjID,'_HRF.pdf']);
print(hrfFigHndl, figNameHrf, '-dpdf', '-r300');