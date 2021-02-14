function [] = analyzeLFContrast_checkResiduals(subjId)
%
% This is to check for systemaic differences in the residuals of the models
% This runs the models and takes the difference in their fits and the
% acutal time course and plots them as a function of stimulus condition
%

display(['STARTING - Main Analysis: ',subjId])
% Load the subject relevant info
analysisParams = getSubjectParams(subjId);


analysisParams.preproc = 'hcp';

analysisParams.saveFigs = false;

% number of TRs to plot for examining the resuduals.
numTR2Plot = 14;

% Number of bootstrap iterations
numIter  = 2;

%turn on or off plotting
analysisParams.showPlots = true;
qcmColor  = [0.4078, 0.2784, 0.5765];
iampColor = [0.8902, 0.6235, 0.5529];

%% Load the relevant data (SDM, HRF, TC)

%set the HRF
[analysisParams] = loadHRF(analysisParams);

% Load the time course
[fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);

% Get a packet for each run (1-20)
[analysisParams, iampTimeCoursePacketPocket,expParams] = generateRunPackets(analysisParams, fullCleanData);

% Concatenate Packet
[analysisParams, theFullPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket,'bootstrap',false);

%% FIT THE TIME COURSE
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% fit the IAMP model
defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
[iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theFullPacket,...
    'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');

% Create the time Course packet
% Get directon/contrast form of time course and IAMP crf packet pockets.
timeCoursePacket = makeDirectionTimeCoursePacketPocket({theFullPacket});

% Make the CRF packet for plotting the stilumi in the nonlinearlty
stimAndRespForPlot = makeDirectionCrfPacketPocket(analysisParams,iampParams);

% Fit the time course with the QCM -- { } is because this expects a cell
[qcmTcOBJ,qcmTcParams] = fitDirectionModel(analysisParams, 'qcmFit', timeCoursePacket,'fitErrorScalar',1000,'talkToMe',false);

% Get the time course predicitions fromt the QCM params fit to the CRF
qcmTimeCourse = responseFromPacket('qcmPred', analysisParams, qcmTcParams{1}, timeCoursePacket, 'plotColor', qcmColor);

%% Get Model Residuals
iampResidual = iampResponses.values - theFullPacket.response.values;
qcmResidual  = qcmTimeCourse{1}.values - theFullPacket.response.values;

%% Scatter Residual
figure; hold on;
scatter(iampResidual,qcmResidual);
axis square;
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
title('Residual Scatter Plot');
ylabel('QCM Residuals');
xlabel('GLM Residuals');
set(gca, 'FontName', 'Helvetica', 'FontSize', 14,'FontWeight', 'normal');

%% Scatter resids with data
figure; hold on;
scatter(theFullPacket.response.values,iampResidual,'MarkerEdgeColor',iampColor);
scatter(theFullPacket.response.values,qcmResidual,'MarkerEdgeColor',qcmColor);
axis square;
hline = refline(-1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
title('Residual Scatter Plot');
ylabel('Residuals');
xlabel('Data');
legend('GLM','QCM')
set(gca, 'FontName', 'Helvetica', 'FontSize', 14,'FontWeight', 'normal');

%% Take the SDM matrix and find the onset times
counter = 0;
meanResidValIAMP= [];
for ii = 1:analysisParams.numSessions
    indx = 1;
    for jj = 1:analysisParams.numAcquisitions
        timingMat = expParams{ii}.timing(:,:,jj);
        
        for kk = 1:size(timingMat,1)
            % The first session here must be the one with the -45,0-,45,90
            if ii == 2
                %shift the coding for the second set of directions
                direction = timingMat(kk,4) + 4;
                contrast  = timingMat(kk,3);
            else
                direction = timingMat(kk,4);
                contrast  = timingMat(kk,3);
            end
            startPt = timingMat(kk,1)+counter;
            stopPt  =timingMat(kk,1)+counter+numTR2Plot;
            if stopPt > length(iampResidual)
                stopPt = length(iampResidual);
            end
            
            residualsDirectionsIAMP{direction}.contrasts(indx,:,contrast) = iampResidual(startPt:stopPt);
            residualsDirectionsQCM{direction}.contrasts(indx,:,contrast)  = qcmResidual(startPt:stopPt);
            
            % get the mean residual for a fixed number of TRs after the
            % stimulus onset (lag)
            lag = 4; %in TRs
            meanResidValIAMP(direction,contrast,indx) = mean(iampResidual(startPt+lag:stopPt));
            meanResidValQCM(direction,contrast,indx) = mean(qcmResidual(startPt+lag:stopPt));
        end
        indx = indx+1;
        counter = counter +360;
    end
end

%% plot the residual plots
for dir = 1:size(analysisParams.directionCoding,2)
    figure
    for con = 1:length(analysisParams.contrastCoding)
        tmpMat = residualsDirectionsIAMP{dir}.contrasts(:,:,con);
        subplot(1,length(analysisParams.contrastCoding),con);
        hold on;
        plot(residualsDirectionsIAMP{dir}.contrasts(:,:,con)','LineWidth',1,'Color',[iampColor,0.6])
        plot(residualsDirectionsQCM{dir}.contrasts(:,:,con)','LineWidth',1,'Color',[qcmColor,0.6])
        plot(mean(residualsDirectionsIAMP{dir}.contrasts(:,:,con),1)','LineWidth',3,'Color',iampColor)
        plot(mean(residualsDirectionsQCM{dir}.contrasts(:,:,con),1)','LineWidth',3,'Color',qcmColor)
        titleStr = getDirContrastString(dir,con);
        title(titleStr);
        ylabel('Residuals');
        xlabel('time (tr)');
        set(gca, 'FontName', 'Helvetica', 'FontSize', 12,'FontWeight', 'normal');
        ylim([-2 2]);
    end
end

%% plot the mean residual as a function of direction (collapsing across contrast)

meanIAMP = mean(mean(meanResidValIAMP,3),2)';
meanQCM  = mean(mean(meanResidValQCM,3),2)';
figure; hold on
angleCoding = [1,5,3,6,2,7,4,8];
scatter(angleCoding,meanIAMP,150,'MarkerEdgeColor',iampColor,...
              'MarkerFaceColor',iampColor,...
              'LineWidth',1.5)
scatter(angleCoding,meanQCM,150,'MarkerEdgeColor',qcmColor,...
              'MarkerFaceColor',qcmColor,...
              'LineWidth',1.5)
ylabel('Mean Residual');
xlabel('Chromatic Dierection');
set(gca,'XTickLabel',{'';'-45';'-22.5';'0';'22.5';'45';'67.5';'90';'112.5';''})
set(gca, 'FontName', 'Helvetica', 'FontSize', 14,'FontWeight', 'normal');
ylim([-.5 .5]);
xlim([0 9]);
hline = refline(0,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('GLM','QCM')




