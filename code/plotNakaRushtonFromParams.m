function [nrVals] = plotNakaRushtonFromParams(amp,exp,semi,varargin)
% Takes in the analysis params struct and returns a voxel by timepoint by
% aquisistion matrix for all the runs found fro the specicied session(s)
%
% Syntax:
%   [fullCleanData, analysisParams, voxelIndex] = getTimeCourse(analysisParams)
%
% Description:
%    This function takes in a struct that is specified in analyzeLFContrast.m
%    and returns a voxel by timepoint by aquisition matrix for all the
%    aquisitions specified in the text files housed in mela_analysis that
%    describe the session(s). f mutliple sessions, they will be
%    concatenated in the 3rd dimension
%
% Inputs:
%    analysisParams    - Stuct contianing relevenat info to the session that
%                        is defined in analyzeLFContrast.m. (string)
%
% Outputs:
%    fullCleanData       - The voxel by timepoint by aquisition matrix
%    analysisParams      - the input analysis params updated with the
%                          number of aquistidiond found per session
%    voxelIndex          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    selectCols
%    plotMovement

% MAB 09/09/18

p = inputParser;
p.addRequired('amp',@isnumeric);
p.addRequired('exp',@isnumeric);
p.addRequired('semi',@isnumeric);
p.addParameter('analysisParams',[],@isstruct)
p.addParameter('xSampleBase',[0:0.01:1],@isnumeric);
p.addParameter('plotFunction',true,@islogical);
p.addParameter('savePlot',false,@islogical);
p.parse(amp,exp,semi,varargin{:})
   
xSampleBase = p.Results.xSampleBase;


NR = @(c) amp .* ((c.^exp)./(c.^exp + semi.^exp));

nrVals = NR(xSampleBase);

if p.Results.plotFunction

    figure
    L1 = plot(xSampleBase,nrVals,'k');
    
    set(L1, 'LineWidth', 2);
    
    hTitle  = title ('Response Nonlinearlity');
    hXLabel = xlabel('Equivalent Contrast'  );
    hYLabel = ylabel('Response');
    
    set(gca,'FontSize',14)
    set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
    set([hXLabel, hYLabel,],'FontSize', 19);
    set( hTitle, 'FontSize', 22,'FontWeight' , 'bold');
    
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'on'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'YTick'       , round(0:max(nrVals)/6:max(nrVals),2) , ...
        'LineWidth'   , 2         , ...
        'ActivePositionProperty', 'OuterPosition');
    
    set(gca, 'Color', 'white' );
end

if p.Results.savePlot
    analysisParams = p.Results.analysisParams;
    figNameCrf =  fullfile(getpref(analysisParams.projectName,'figureSavePath'),analysisParams.expSubjID, ...
    [analysisParams.expSubjID,'_Nonlinearity_' analysisParams.sessionNickname '_' analysisParams.preproc '.pdf']);
    FigureSave(figNameCrf,gca,'pdf');
end
