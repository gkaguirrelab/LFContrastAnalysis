function LFContrastAnalysisLocalHook
%  LFContrastAnalsysisLocalHook
%
% Configure things for working on the  LFContrastAnalysis project.
%
% For use with the ToolboxToolbox.
%
% If you 'git clone' ILFContrastAnalsysis into your ToolboxToolbox "projectRoot"
% folder, then run in MATLAB
%   tbUseProject('LFContrastAnalysis')
% ToolboxToolbox will set up IBIOColorDetect and its dependencies on
% your machine.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/LFContrastAnalsysisLocalHook.m
%
% Each time you run tbUseProject('LFContrastAnalysis'), ToolboxToolbox will
% execute your local copy of this file to do setup for LFContrastAnalsysis.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


%% Say hello.
fprintf('LFContrastAnalysis local hook.\n');
projectName = 'LFContrastAnalysis';

%% Delete any old prefs
if (ispref(projectName))
    rmpref(projectName);
end

%% Specify base paths for materials and data
[~, userID] = system('whoami');
userID = strtrim(userID);
switch userID
    case {'dhb'}
        melaMaterialsPath = ['/Users1/DropboxLab/MELA_materials'];
        melaDatabasePath  = ['/Users1/DropboxLab/MELA_data/'];
        melaAnalysisPath  = ['/Users1/DropboxLab/MELA_analysis/'];
        figureSavePath    = fullfile(melaAnalysisPath,projectName,'Figures');
    case {'michael'}
        melaMaterialsPath = ['/Users/' userID '/labDropbox/MELA_materials'];
        melaDatabasePath  = ['/Users/' userID '/labDropbox/MELA_data/'];
        melaAnalysisPath  = ['/Users/' userID '/labDropbox/MELA_analysis/'];
        figureSavePath    = fullfile(melaAnalysisPath,projectName,'Figures');
    otherwise
        melaMaterialsPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        melaDatabasePath  = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
        melaAnalysisPath  = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/'];
        figureSavePath    = fullfile(melaAnalysisPath,projectName,'Figures');
end

%% Specify where output goes

if ismac
    % Code to run on Mac plaform
    setpref(projectName,'analysisScratchDir','/tmp/flywheel');
    setpref(projectName,'projectRootDir',fullfile('/Users/',userID,'/Documents/flywheel',projectName));
    setpref(projectName,'projectPath', fullfile(melaDatabasePath,'Experiments','OLApproach_TrialSequenceMR'));
    setpref(projectName,'melaAnalysisPath', melaAnalysisPath);
    setpref(projectName,'figureSavePath', figureSavePath);
    setpref(projectName,'materialsPath',fullfile(melaMaterialsPath,'Experiments','OLApproach_TrialSequenceMR'));
elseif isunix
    % Code to run on Linux plaform
    setpref(projectName,'analysisScratchDir','/tmp/flywheel');
    setpref(projectName,'projectRootDir',fullfile('/home/',userID,'/Documents/flywheel',projectName));
    setpref(projectName,'projectPath', fullfile(melaDatabasePath,'Experiments','OLApproach_TrialSequenceMR'));
    setpref(projectName,'melaAnalysisPath', melaAnalysisPath);
    setpref(projectName,'figureSavePath', figureSavePath);
    setpref(projectName,'materialsPath',fullfile(melaMaterialsPath,'Experiments','OLApproach_TrialSequenceMR'));
elseif ispc
    % Code to run on Windows platform
    warning('No supported for PC')
else
    disp('What are you using?')
end
