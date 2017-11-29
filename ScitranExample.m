


%% Open a scitran object
%
% Maybe this worked
st = scitran('upennflywheel');

%% Do we have access to any projects
projects = st.search('projects');

%% Find sessions in the first project
projectID = projects{1}.project.x_id;
sessions = st.search('session',...
    'project id',projectID,...
    'summary',true);

analyses = st.search('analysis',...
    'session label contains','2017-09-19 14:35',...
    'summary',true);

%% Try to download something
outfile = st.downloadObject('5a15c112e108ff001bc0d2c1');

%% Upload a file
fullFileName = fullfile('/Users/dhb/Desktop/ByeByeSensorNotes.txt');
st.upload(fullFileName,'session',sessions{1}.session.x_id);
files = st.search('file',...
      'file name','ByeByeSensorNotes.txt');


%% Specify what we're after
project       = 'LFContrast';
session_label = '2017-09-19 14:35';
subject_code  = 'HERO_gka1';

%% Some client init code that we don't understand
%st.toolbox('toolboxes.json','project',project);

%% Get a file the scitran way
theFiles = st.search('files',...
                               'project label',project,...
                               'session label',session_label);
% fnameElectrodes = fullfile(workDir,'sub-19_loc.tsv');
% st.get(electrodePositions{1},'destination',fnameElectrodes);


%% Put a file the scitran way
%
% Query for the target
target = st.search('sessions',...
                   'project label', project,...
                   'subject code', subject_code);
% Build uplaod struct with the information needed to upload the result
upload.file = fullfile(workDir,[subject_code, '_ElectrodePositionsAndLabels.png']);
upload.containerType = 'sessions';
upload.id = target{1}.id;

% Do the upload
st.put('file', upload);