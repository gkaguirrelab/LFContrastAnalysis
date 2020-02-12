% This is a script to determine the difference in fitting the models to the
% whole time course versus the average run by run fits.

% Initialize
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
subjId = 'KAS25';

% Fit the model to the concatenated time series 
[modelResponseStructIAMP, modelResponseStructQCM, thePacketIAMP, thePacketQCM] = analyzeLFContrast_fullTimeSeries(subjId,'showPlot',false);

% Calculate R^2
corrValsIAMP = [modelResponseStructIAMP.values',thePacketIAMP.response.values'];
rSquaredIAMP = corr(corrValsIAMP).^2;

% Calculate R^2
corrValsQCM = [modelResponseStructQCM.values',thePacketIAMP.response.values'];
rSquaredQCM = corr(corrValsQCM).^2;


% CHOP UP TIME COURSE PREDICTIONS TO MATCH THE RUN LENGTH

%% Fitting with average run by run method




















% plot it


