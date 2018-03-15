function []= plotTimeCourse(timeCourse,block,responseStruct)
% using this to look at stuff not meant for anything more right now...
% i will turn this into a proper function with comments soon 



%% Plotting the power levels and trial starts, trial stops, and trial wait times.
for i = 1:length(responseStruct.events)
    % set up trial start and stop time markers
    trialStartTime(i) = responseStruct.events(i).tTrialStart - responseStruct.tBlockStart;
    trialEndTime(i) = responseStruct.events(i).tTrialEnd - responseStruct.tBlockStart;
    trialWaitTime(i) = trialStartTime(i) + responseStruct.events(i).trialWaitTime;
    
    % set up power level plot relevent vars.
    timeStep = block(i).modulationData.modulationParams.timeStep;
    stimulusDuration = block(i).modulationData.modulationParams.stimulusDuration;
    sampleBasePowerLevel{i} =  (trialStartTime(i) + responseStruct.events(i).trialWaitTime):timeStep:(trialStartTime(i) + responseStruct.events(i).trialWaitTime+ stimulusDuration -timeStep);
end

figure;
subplot(3,1,1); hold on;
title('Power Level Modulations') 
for ii = 1:length(responseStruct.events)
    plot([trialStartTime(ii) trialStartTime(ii)]-0.1, [-1 1],'r','LineWidth',2);
    plot([trialEndTime(ii) trialEndTime(ii)], [-1 1],'b','LineWidth',2);
    plot([trialWaitTime(ii) trialWaitTime(ii) ], [-1 1],'g--');
    plot(sampleBasePowerLevel{ii},block(ii).modulationData.modulation.powerLevels,'k');
end


subplot(3,1,2); hold on;
title('Time Course') 
timepoints = 0.8.*[1:size(timeCourse)]-0.8;
plot(timepoints,timeCourse);
matVals = timeCourse(:,3:end);
minVal = round(min(matVals(:)) - 100);
maxVal = round(max(matVals(:)) + 100);
for ii = 1:length(responseStruct.events)
    plot([trialStartTime(ii) trialStartTime(ii)]-0.1, [minVal maxVal],'r','LineWidth',2);
    plot([trialEndTime(ii) trialEndTime(ii)], [minVal maxVal],'b','LineWidth',2);
    plot([trialWaitTime(ii) trialWaitTime(ii) ], [minVal maxVal],'g--');
end

subplot(3,1,3); hold on;
title('percent signal change (relative to mean)') 
timepoints = 0.8.*[1:size(timeCourse)]-0.8;
meanMat = repmat(mean(timeCourse,1),[size(timeCourse,1),1]);

PSC = 100.*((timeCourse - meanMat)./meanMat);
plot(timepoints,PSC);
y = plot(timepoints,mean(PSC,2),'--k');
for ii = 1:length(responseStruct.events)
    plot([trialStartTime(ii) trialStartTime(ii)]-0.1, [-10 10],'r','LineWidth',2);
    plot([trialEndTime(ii) trialEndTime(ii)], [-10 10],'b','LineWidth',2);
    plot([trialWaitTime(ii) trialWaitTime(ii) ], [-10 10],'g--');
end
ylim([-3 3])
yticks([-3 -2 -1 0 1 2 3])
set(gca, 'YGrid', 'on', 'XGrid', 'off')

set(y,'LineWidth',2)
legend(y,'Mean PSC')
end