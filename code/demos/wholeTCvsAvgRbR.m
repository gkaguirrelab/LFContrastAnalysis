% This is a script to determine the difference in fitting the models to the
% whole time course versus the average run by run fits.

% Initialize
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
subjId = 'KAS25';

% Fit the model to the concatenated time series
[modelResponseStructIAMP_FTC, modelResponseStructQCM_FTC, thePacketIAMP_FTC, thePacketQCM_FTC] = fullTimeCourseFit(subjId);

[modelResponseStructIAMP_RbR, modelResponseStructQCM_RbR, thePacketIAMP_RbR, thePacketQCM_RbR]  = analyzeLFContrast_runByRun(subjId);

QCM_FTC_ChoppedUp = chopUpTimeCourse(modelResponseStructQCM_FTC,20);
timeCourses = [thePacketIAMP_RbR(1,1:end),thePacketIAMP_RbR(2,1:end)];


for ii = 1:2
    figure;hold on
    for jj = 1:10
        indx = ((ii-1)*10) + jj;
        
        % Calculate R^2
        corrValsQCM_FTC = [QCM_FTC_ChoppedUp{indx}.values',timeCourses{indx}.response.values'];
        rSquaredQCM_FTC = corr(corrValsQCM_FTC).^2;
        rSquaredQCM_FTC_vals(indx) = rSquaredQCM_FTC(2);
        % Calculate R^2
        corrValsQCM_RbR = [modelResponseStructQCM_RbR{indx}.values',timeCourses{indx}.response.values'];
        rSquaredQCM_RbR = corr(corrValsQCM_RbR).^2;
        rSquaredQCM_RbR_vals(indx) = rSquaredQCM_RbR(2);
        
        % plot it
        subplot(4,3,jj)
        hold on
        plot(QCM_FTC_ChoppedUp{indx}.timebase,QCM_FTC_ChoppedUp{indx}.values,'r')
        plot(modelResponseStructQCM_RbR{indx}.timebase,modelResponseStructQCM_RbR{indx}.values,'b')
        plot(timeCourses{indx}.response.timebase,timeCourses{indx}.response.values,'k')
        title(sprintf('run %s: r2-ftc =%s, r2-rbr =%s',num2str(indx),num2str(rSquaredQCM_FTC(2)),num2str(rSquaredQCM_RbR(2))))
    end
    legend('Full Time Course','Median of Run Fits', 'Time Course')
end
display(sprintf(' Average R^2 for FTC = %s',num2str(mean(rSquaredQCM_FTC_vals))))
display(sprintf(' Average R^2 for RbR = %s',num2str(mean(rSquaredQCM_RbR_vals))))
