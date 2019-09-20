% Initialize
clear;

%% Set up bootstrap
%
% number of bootstrap iterations
numIterations = 10;
% Subject: 'LZ23', 'KAS25', 'AP26'
subj = 'AP26';
percentile = 65;

%% Session 1 Bootstrap
%
% Get subject specific params
analysisParams = getSubjectParams(subj);

% Get the time course data for the subject
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Fit the GLM 
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% Bootstrap the params
[boot_qcmParamsOrig,boot_nrValsOrig] = fit_bootstrapQCM(analysisParams,iampParams,numIterations);

% get the median and 65th percentile
sortedRows = sortrows(boot_nrValsOrig,size(boot_nrValsOrig,2));

if mod(size(sortedRows,1),2) == 0 
    medianNROrig = mean([sortedRows(size(sortedRows,1)/2,:);sortedRows(1+size(sortedRows,1)/2,:)]);
else
    medianNROrig = sortedRows(ceil(size(sortedRows,1)/2),:);
end

% get error bars
errorIndx = (numIterations-((percentile/100)*numIterations))/2;
if floor(errorIndx) == errorIndx
    nrOrigUB = sortedRows(end-errorIndx,:);
    nrOrigLB = sortedRows(errorIndx,:);
    nrOrigErrorBars = [nrOrigUB-medianNROrig;medianNROrig-nrOrigLB];
else
    nrOrigUB = mean([sortedRows(end-ceil(errorIndx),:);sortedRows(end-floor(errorIndx),:)]);
    nrOrigLB = mean([sortedRows(ceil(errorIndx),:);sortedRows(floor(errorIndx),:)]);
    nrOrigErrorBars = [nrOrigUB-medianNROrig;medianNROrig-nrOrigLB];
end

%% Get the model fit 
[qcmParamsOrig, nrValsOrig] = fit_QCM(analysisParams,iampParams);

%% Replication Bootstrap
%
% Get subject specific params
analysisParams = getSubjectParams([subj '_replication']);

% Get the time course data for the subject
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Fit the GLM 
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% Bootstrap the params
[boot_qcmParamsRep,boot_nrValsRep] = fit_bootstrapQCM(analysisParams,iampParams,numIterations);

% get the median and 65th percentile
sortedRows = sortrows(boot_nrValsRep,size(boot_nrValsRep,2));

if mod(size(sortedRows,1),2) == 0 
    medianNRRep = mean([sortedRows(size(sortedRows,1)/2,:);sortedRows(1+size(sortedRows,1)/2,:)]);
else
    medianNRRep = sortedRows(ceil(size(sortedRows,1)/2),:);
end

% get error bars
errorIndx = (numIterations-((percentile/100)*numIterations))/2;
if floor(errorIndx) == errorIndx
    nrRepUB = sortedRows(end-errorIndx,:);
    nrRepLB = sortedRows(errorIndx,:);
    nrRepErrorBars = [nrRepUB-medianNRRep;medianNRRep-nrRepLB];
else
    nrRepUB = mean([sortedRows(end-ceil(errorIndx),:);sortedRows(end-floor(errorIndx),:)]);
    nrRepLB = mean([sortedRows(ceil(errorIndx),:);sortedRows(floor(errorIndx),:)]);
    nrRepErrorBars = [nrRepUB-medianNRRep;medianNRRep-nrRepLB];
end


%% Get the model fit 
[qcmParamsRep, nrValsRep] = fit_QCM(analysisParams,iampParams)


%% Plotting
figure; hold on
scatter(boot_qcmParamsOrig(1,:),boot_qcmParamsOrig(2,:),'r')
scatter(boot_qcmParamsRep(1,:),boot_qcmParamsRep(2,:),'b')
hTitle  = title ('Parameter Bootstrap Scatter Plot');
hXLabel = xlabel('Minor Axis Ratio'  );
hYLabel = ylabel('Angle (degrees)');
set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 14);
set( hTitle, 'FontSize', 14,'FontWeight' , 'bold');

figure; hold on
shadedErrorBars(0:0.01:1,medianNROrig,nrOrigErrorBars,'lineProps','r')
shadedErrorBars(0:0.01:1,medianNRRep,nrRepErrorBars,'lineProps','b')
plot(0:0.01:1,nrValsOrig,'--r','linewidth',2)
plot(0:0.01:1,nrValsRep,'--b','linewidth',2)

hTitle  = title ('Nonlinearity Bootstrap');
hXLabel = xlabel('Eq. Contrast'  );
hYLabel = ylabel('Response');


set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 14);
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
    'YTick'       , 0:2:10 , ...
    'LineWidth'   , 2         , ...
    'ActivePositionProperty', 'OuterPosition');
ylim([0 10]);

set(gcf, 'Color', 'white' );