%% setupLFContrastData
%
% Script that downloads the analysis products from flywheel into the
% directory specified in the local hook file and organizes it in a
% standard format.

% MAB 03/2018

%% Convenience variables
projectName = 'LFContrastAnalysis';
flywheelName = 'LFContrast';

%% Analysis labels that we are going to go and get
neuropythyLabel = 'retinotopy-templates 10/13/2018 19:28:00'

%% Make Benson Gear (Neuropythy) directory
bensonDir = fullfile('/Users/michael/Downloads/neuropythy');
if ~exist(bensonDir)
    mkdir(bensonDir)
end

%% Download Benson Gear (Neuropythy) output
clear fwInfo
[fwInfo] = getAnalysisFromFlywheel(flywheelName,neuropythyLabel,bensonDir, 'verbose', true);
