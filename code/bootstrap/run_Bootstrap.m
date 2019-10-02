% Initialize
clear;

%% Set up bootstrap
%
% number of bootstrap iterations
numIterations = 100;
% Subject: 'LZ23', 'KAS25', 'AP26'
subj = 'LZ23';
percentile = 65;

%  #############################################
%% ############ Session 1  Bootstrap ###########
%  #############################################

% Get subject specific params
analysisParams = getSubjectParams(subj);

% Get the time course data for the subject
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Fit the GLM
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

% Bootstrap the params
[boot_qcmParamsOrig,boot_nrValsOrig] = fit_bootstrapQCM(analysisParams,iampParams,numIterations);

% get the median
sortedRows = sortrows(boot_nrValsOrig,size(boot_nrValsOrig,2));
if mod(size(sortedRows,1),2) == 0
    medianNROrig = mean([sortedRows(size(sortedRows,1)/2,:);sortedRows(1+size(sortedRows,1)/2,:)]);
else
    medianNROrig = sortedRows(ceil(size(sortedRows,1)/2),:);
end

% get the CI
errorIndx = (numIterations-((percentile/100)*numIterations))/2;
for ff = 1:size(boot_nrValsOrig,2)
    sortedTimePointVals = sort(boot_nrValsOrig(:,ff));
    if floor(errorIndx) == errorIndx
        nrOrigUB(ff) = sortedTimePointVals(end-errorIndx,:);
        nrOrigLB(ff) = sortedTimePointVals(errorIndx,:);
    else
        nrOrigUB(ff) = mean([sortedTimePointVals(end-ceil(errorIndx),:);sortedTimePointVals(end-floor(errorIndx),:)]);
        nrOrigLB(ff) = mean([sortedTimePointVals(ceil(errorIndx),:);sortedTimePointVals(floor(errorIndx),:)]);
    end
end
nrOrigErrorBars = [nrOrigUB-medianNROrig;medianNROrig-nrOrigLB];

%% Get the model fit
[qcmParamsOrig, nrValsOrig] = fit_QCM(analysisParams,iampParams);

%% Get Median for ellipse params
sortedMinorAxis = sort(boot_qcmParamsOrig(1,:))';
sortedAngle     = sort(boot_qcmParamsOrig(2,:))';
errorIndx = (length(sortedAngle)-((percentile/100)*length(sortedAngle)))/2;
if mod(size(sortedMinorAxis,1),2) == 0
    medianAngleOrig= mean([sortedAngle(size(sortedAngle,1)/2,:);sortedAngle(1+size(sortedAngle,1)/2,:)]);
    medianMarOrig = mean([sortedMinorAxis(size(sortedMinorAxis,1)/2,:);sortedMinorAxis(1+size(sortedMinorAxis,1)/2,:)]);
else
    medianAngleOrig = sortedAngle(ceil(size(sortedAngle,1)/2),:);
    medianMarOrig = sortedMinorAxis(ceil(size(sortedMinorAxis,1)/2),:);
end
%% Get CI for ellipse params
if floor(errorIndx) == errorIndx
    angleRepUBorig = sortedAngle(end-errorIndx,:);
    angleRepLBorig = sortedAnlge(errorIndx,:);
    marRepUBorig = sortedMinorAxis(end-errorIndx,:);
    marRepLBorig = sortedMinorAxis(errorIndx,:);
else
    angleRepUBorig = mean([sortedAngle(end-ceil(errorIndx),:);sortedAngle(end-floor(errorIndx),:)]);
    angleRepLBorig = mean([sortedAngle(ceil(errorIndx),:);sortedAngle(floor(errorIndx),:)]);
    marRepUBorig = mean([sortedMinorAxis(end-ceil(errorIndx),:);sortedMinorAxis(end-floor(errorIndx),:)]);
    marRepLBorig = mean([sortedMinorAxis(ceil(errorIndx),:);sortedMinorAxis(floor(errorIndx),:)]);
end

%  #############################################
%% ############ Replication Bootstrap ##########
%  #############################################

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
for ff = 1:size(boot_nrValsRep,2)
    sortedTimePointVals = sort(boot_nrValsRep(:,ff));
    if floor(errorIndx) == errorIndx
        nrRepUB(ff) = sortedTimePointVals(end-errorIndx,:);
        nrRepLB(ff) = sortedTimePointVals(errorIndx,:);
    else
        nrRepUB(ff) = mean([sortedTimePointVals(end-ceil(errorIndx),:);sortedTimePointVals(end-floor(errorIndx),:)]);
        nrRepLB(ff) = mean([sortedTimePointVals(ceil(errorIndx),:);sortedTimePointVals(floor(errorIndx),:)]);
        
    end
end
nrRepErrorBars = [nrRepUB-medianNRRep;medianNRRep-nrRepLB];

%% Get the model fit
[qcmParamsRep, nrValsRep] = fit_QCM(analysisParams,iampParams)

%% Get Median for ellipse params
sortedMinorAxis = sort(boot_qcmParamsRep(1,:))';
sortedAngle     = sort(boot_qcmParamsRep(2,:))';
errorIndx = (length(sortedAngle)-((percentile/100)*length(sortedAngle)))/2;
if mod(size(sortedMinorAxis,1),2) == 0
    medianAngleRep = mean([sortedAngle(size(sortedAngle,1)/2,:);sortedAngle(1+size(sortedAngle,1)/2,:)]);
    medianMarRep= mean([sortedMinorAxis(size(sortedMinorAxis,1)/2,:);sortedMinorAxis(1+size(sortedMinorAxis,1)/2,:)]);
else
    medianAngleRep = sortedAngle(ceil(size(sortedAngle,1)/2),:);
    medianMarRep= sortedMinorAxis(ceil(size(sortedMinorAxis,1)/2),:);
end
%% Get CI for ellipse params
if floor(errorIndx) == errorIndx
    angleRepUBrep = sortedAngle(end-errorIndx,:);
    angleRepLBrep = sortedAnlge(errorIndx,:);
    marRepUBrep = sortedMinorAxis(end-errorIndx,:);
    marRepLBrep = sortedMinorAxis(errorIndx,:);
else
    angleRepUBrep = mean([sortedAngle(end-ceil(errorIndx),:);sortedAngle(end-floor(errorIndx),:)]);
    angleRepLBrep = mean([sortedAngle(ceil(errorIndx),:);sortedAngle(floor(errorIndx),:)]);
    marRepUBrep = mean([sortedMinorAxis(end-ceil(errorIndx),:);sortedMinorAxis(end-floor(errorIndx),:)]);
    marRepLBrep = mean([sortedMinorAxis(ceil(errorIndx),:);sortedMinorAxis(floor(errorIndx),:)]);
end


%% Plotting
figure; hold on
scatter(boot_qcmParamsOrig(1,:),boot_qcmParamsOrig(2,:),'r')
errorbar(medianMarOrig,medianAngleOrig,abs(angleRepLBorig-medianAngleOrig),...
         abs(angleRepUBorig-medianAngleOrig),abs(medianMarOrig-marRepLBorig),abs(medianMarOrig-marRepUBorig),'*r')
     
scatter(boot_qcmParamsRep(1,:),boot_qcmParamsRep(2,:),'b')
errorbar(medianMarRep,medianAngleRep,abs(angleRepLBrep-medianAngleRep),...
         abs(angleRepUBrep-medianAngleRep),abs(medianMarRep-marRepLBrep),abs(medianMarRep-marRepUBrep),'*b')
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