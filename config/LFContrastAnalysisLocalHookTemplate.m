function LFContrastAnalsysisLocalHook
%  LFContrastAnalsysisLocalHook
%
% Configure things for working on the  LFContrastAnalsysis project.
%
% For use with the ToolboxToolbox.
%
% If you 'git clone' ILFContrastAnalsysis into your ToolboxToolbox "projectRoot"
% folder, then run in MATLAB
%   tbUseProject('LFContrastAnalsysis')
% ToolboxToolbox will set up IBIOColorDetect and its dependencies on
% your machine.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/LFContrastAnalsysisLocalHook.m
%
% Each time you run tbUseProject('LFContrastAnalsysis'), ToolboxToolbox will
% execute your local copy of this file to do setup for LFContrastAnalsysis.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


%% Say hello.
fprintf('LFContrastAnalsysis local hook.\n');
projectName = 'LFContrastAnalsysis';

% %% UnitTestToolbox and RemoteDataToolbox setup.
% %
% % If you customize your rdt-config json file, you will want to place it
% % somewhere outside the repository and change the path to point to your
% % copy.  The usual reason for doing this is so that you can add your
% % username and password and have write permission.
% projectBaseDir = tbLocateProject(projectName);
% rdtConfig = fullfile(projectBaseDir, 'configuration', ['rdt-config-' projectName '.json']);
% p = struct(...
%     'projectName',           projectName, ...                                                                                 % The project's name (also the preferences group name)
%     'validationRootDir',     IBIOCDValidationDir, ...                                                                         % Directory location where the 'scripts' subdirectory resides.
%     'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
%     'alternateFullDataDir',  '', ...                                                                                          % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
%     'useRemoteDataToolbox',  true, ...                                                                                        % If true use Remote Data Toolbox to fetch full validation data on demand.
%     'remoteDataToolboxConfig', rdtConfig, ...                                                                                 % Struct, file path, or project name with Remote Data Toolbox configuration.
%     'clonedWikiLocation',    '', ...                                                                                          % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
%     'clonedGhPagesLocation', '', ...                                                                                          % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
%     'githubRepoURL',         '', ...                                                                                          % Github URL for the project. This is only used for publishing tutorials.
%     'generateGroundTruthDataIfNotFound',true,...                                                                              % Flag indicating whether to generate ground truth if one is not found
%     'listingScript',         'IBIOCDValidateListAllValidationDirs', ...                                                       % Script that lists dirs to find validation scripts in
%     'coreListingScript',     '', ...                                                                                          % Not used in this project
%     'numericTolerance',      1e-11 ...                                                                                        % Numeric tolerance for comparisons with validation data.
%     );
% 
% generatePreferenceGroup(p);
% UnitTest.usePreferencesForProject(p.projectName);
% 
% 
% %% Output directory.
% %
% % This is where the project writes its output.  By default, we'll stick it
% % in a subfolder of a folder called output, in the tbUserFolder.  But you
% % may want it somewhere else.
% outputBaseDir = fullfile(tbUserFolder(), 'output', projectName);
% if (7 ~= exist(outputBaseDir, 'dir'))
%     mkdir(outputBaseDir);
% end
% 
% setpref(projectName, 'outputBaseDir', outputBaseDir);
% 
% end
% 
% 
% %% Generate preferences that work with UnitTest toolbox.
% function generatePreferenceGroup(p)
% % Remove any existing preferences for this project
% if ispref(p.projectName)
%     rmpref(p.projectName);
% end
% 
% % Renerate and save the project-specific preferences
% setpref(p.projectName, 'projectSpecificPreferences', p);
% fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
% end
