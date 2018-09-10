 % Convenience variables
analysisParams.projectName       = 'LFContrastAnalysis';
analysisParams.flywheelName      = 'LFContrast';
analysisParams.subjID            = 'sub-HEROGKA1';
analysisParams.expSubjID         = 'HERO-gka1';
analysisParams.session           = 'ses-ResearchAguirre';
analysisParams.sessionFolderName = 'HERO_GKA1_2018-07-28';
analysisParams.sessionDate       = '2018-07-28';
analysisParams.sessionNumber     = 'session_1';
analysisParams.sessionDir = fullfile(getpref('LFContrastAnalysis','projectRootDir'),sessionFolderName);


% make mask from the area and eccentricity maps
analysisParams.areaNum     = 2;
analysisParams.eccenRange  = [3 20];

% Define the TR
analysisParams.TR = 0.800;

getTimeCourse(analysisParams)


%% Get trial order info:
trialOrderDir = fullfile(getpref(analysisParams.projectName,'melaDataPath'),HERO_gka1/2018-07-28/session_1');
trialOrderFiles = {'CRF_analysisParams.session_1_scan1.mat' ...
    'CRF_session_1_scan2.mat', ...
    'CRF_session_1_scan3.mat', ...
    'CRF_session_1_scan4.mat', ...
    'CRF_session_1_scan5.mat', ...
    'CRF_session_1_scan6.mat', ...
    'CRF_session_1_scan7.mat', ...
    'CRF_session_1_scan8.mat', ...
    'CRF_session_1_scan9.mat', ...
    'CRF_session_1_scan10.mat'};