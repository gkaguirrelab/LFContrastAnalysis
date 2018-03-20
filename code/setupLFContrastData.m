%% setupLFContrastData
%
% MAB: Write a script here to download the analysis products from flywheel
% into the  directory specified in the local hook file and organize it in a
% nice way

%% Check for project level directory 
projectDir = getpref('LFContrastAnalysis','projectRootDir');
if ~exist(projectDir)
    mkdir(projectDir)
end

%% Set up fmriprep directory and download

fmriprepDir = fullfile(projectDir,'fmriprep')

if ~exist(fmriprepDir)
    mkdir(fmriprepDir)
end

% download the data 
theProject    = 'LFContrast';
analysisLabel = 'fmriprep 02/09/2018 11:40:55';
analysis_id = getAnalysisFromFlywheel(theProject,analysisLabel,fmriprepDir, 'verbose', true)

