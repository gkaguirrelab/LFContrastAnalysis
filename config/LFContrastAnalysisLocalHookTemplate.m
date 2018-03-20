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
    case {'melanopsin' 'pupillab'}
        materialsBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        MELA_dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
    case {'dhb'}
        materialsBasePath = ['/Users1' '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        MELA_dataBasePath = ['/Users1' '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];     
    case {'nicolas'}
        materialsBasePath = '/Volumes/Manta TM HD/Dropbox (Aguirre-Brainard Lab)/MELA_materials';
        MELA_dataBasePath = '/Volumes/Manta TM HD/Dropbox (Aguirre-Brainard Lab)/MELA_data';
    otherwise
        materialsBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_materials'];
        MELA_dataBasePath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
end

%% Specify where output goes
setpref(projectName,'analysisScratchDir','/tmp/flywheel');
setpref(projectName,'projectRootDir',fullfile('/Users/',userID,'/Documents/flywheel',projectName));

