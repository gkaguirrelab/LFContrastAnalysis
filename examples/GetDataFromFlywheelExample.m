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
    sessionAcqs{ii} = fw.getSessionAcquisitions(sessionId{ii})
    sessionAnalyses{ii} = fw.getSessionAnalysis(sessionId{ii})
end

