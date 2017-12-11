%% Analyze LFContrast Data.
%
% This script calls function in order to analyze the data for the
% LFContrast experiment. 

%% Set up params.
params.subjID       = 'HEROgka1';
params.session      = 'ses-201709191435';
params.regFileName  = 'register.dat';
params.atlasInName  = 'inMap.nii.gz'
params.atlasOutName = 'outMap.nii.gz'

%% Align Atlas to Func
path2ResampAtlas = applyRetAtlas2Functional(retFiles,path2FuncData,params)

%% Extract Signal from voxels
meanSignal = extractMeanSignalFromROI(timeSeries,ROI)

%% Plot the data
figure; hold on
