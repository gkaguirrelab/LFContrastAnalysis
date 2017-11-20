% GetDataFromFlywheelExample  Use api to snag some data from flywheel
%
% Description:
%   Use our machinery to snag a sample dataset from flywheel
%

% History
%  11/10/17  dhb, gka, mab  Wrote it because we are so excited.

%% Clear
clear; close all;

%% Define our project
theProject = 'LFContrast';

%% Open flywheel object
fw = Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

%% Find out who we are
me = fw.getCurrentUser();
fprintf('I am %s %s\n', me.firstname, me.lastname);

%% Get a list of our projects
theProjectIndex = [];
projects = fw.getAllProjects();
%fprintf('Avaliable projects\n');
for ii = 1:length(projects)
    %fprintf('\t%s\n',projects{ii}.label)
    if (strcmp(theProject,projects{ii}.label))
        theProjectIndex = ii;
        break;
    end
end
if (isempty(theProjectIndex))
    error('Could not find specified project %s\n',theProject);
end
fprintf('Found project %s!\n',projects{theProjectIndex}.label);
projectId = projects{ii}.id;

%% Try to get output from fmriPrep for each session
projectSessions = fw.getProjectSessions(projectId);
for ii = 1:length(projectSessions)
    % Get session ID
    if (length(projectSessions) == 1)
        sessionId{ii} = projectSessions.id;
    else
        sessionId{ii} = projectSessions{ii}.id;
    end
    
    % Get acquisitions for each session
    sessionAcqs{ii} = fw.getSessionAcquisitions(sessionId{ii});
end

%% Try to download the output of an analysis
%
% Given some analysis label, download the files that were generated.
%
% For this we use: 
%   fw.downloadFileFromAnalysis(session_id, analysis_id, file_name, output_name)

% The label for the analysis 
analysis_label = 'fmriprep 10/26/2017 22:17:09';

% Where do you want the files stored? 
out_dir = getpref('LFContrastAnalysis','analysisScratchDir');
if (~exist(out_dir,'dir'))
    mkdir(out_dir);
end

%% Set-up search structure and search 
searchStruct = struct('return_type', 'file', ...
    'filters', {{struct('term', ...
    struct('analysis0x2elabel', analysis_label))}});
results = fw.search(searchStruct);

methods Flywheel
fw.getSdkVersion

% Iterate over results and download the files
for ii = 1:numel(results)
    file_name = results(ii).file.name;
    output_name = fullfile(out_dir, file_name);
    
    session_id = results(ii).session.x_id;
    analysis_id = results(ii).analysis.x_id;

    fprintf('Downloading %dMB file: %s ... \n', round(results(ii).file.size / 1000000), file_name);
    tic; fw.downloadFileFromAnalysis(session_id, analysis_id, file_name, output_name); toc
end
