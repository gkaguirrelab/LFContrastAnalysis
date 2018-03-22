%% setupLFContrastData
%
% Script that downloads the analysis products from flywheel into the 
% directory specified in the local hook file and organizes it in a
% standard format.

% MAB 03/2018

%% Check for/make the project level directory 
projectDir = getpref('LFContrastAnalysis','projectRootDir');
if ~exist(projectDir)
    mkdir(projectDir)
end

%% Set up fmriprep directory and download

fmriprepDir = fullfile(projectDir,'fmriprep');

if ~exist(fmriprepDir)
    mkdir(fmriprepDir)
end

% download the data 
theProject    = 'LFContrast';
analysisLabel = 'fmriprep 02/09/2018 11:40:55';
[fwInfo] = getAnalysisFromFlywheel(theProject,analysisLabel,fmriprepDir, 'verbose', true, 'searchDir', projectDir);


%% make session dir
sessionDate = fwInfo.timestamp(1:10);
sessionFolder = [fwInfo.subject,'_',sessionDate];
sessionDir = fullfile(projectDir,sessionFolder);

if ~exist(sessionDir)
    mkdir(sessionDir)
end

%% Move contents of fmriprep into session directory 
dataLocation = fullfile(fmriprepDir, [fwInfo.subject, '_', fwInfo.analysis_id],fwInfo.analysis_id);
% move fmriprep output
if exist(fullfile(dataLocation,'fmriprep'))
    movefile(fullfile(dataLocation,'fmriprep'), sessionDir)
else
    warning('LFContrastAnalysis:analysisDirectoryMissing','WARNING: fmriprep directory missing from %s folder',fwInfo.analysis_id);
end

% move freesurfer output
if exist(fullfile(dataLocation,'freesurfer'))
   movefile(fullfile(dataLocation,'freesurfer'), sessionDir)
else
    warning('LFContrastAnalysis:analysisDirectoryMissing','WARNING: freesurfer directory missing from %s folder',fwInfo.analysis_id);
end

%% Make bookkeping files

% create bookkeeping file so that the download function can find the
% analysis id after the clean up at the end of this script. This allows 
% for download skip

% fmriprep file
cmd = ['echo -e "' 'This is a bookkeeping file.\nDo not delete unless you delete the whole analysis directory.\nIf you want to download the data again and you got both frmiprep and freesurfer from the same gear, then you need to dete both analysis directories" > ' ...
       fullfile(sessionDir,'fmriprep',[fwInfo.analysis_id '.txt'])];
system(cmd);

% freesurfer file
cmd = ['echo -e "' 'This is a bookkeeping file.\nDo not delete unless you delete the whole analysis directory.\nIf you want to download the data again and you got both frmiprep and freesurfer from the same gear, then you need to dete both analysis directories" > ' ...
       fullfile(sessionDir,'freesurfer',[fwInfo.analysis_id '.txt'])];
system(cmd);


%% Make Benson Gear (Neuropythy) directory
bensonDir = fullfile(sessionDir,'neuropythy');
if ~exist(bensonDir)
    mkdir(bensonDir)
end

%% Download Benson Gear (Neuropythy) output
clear fwInfo
theProject    = 'LFContrast';
analysisLabel = 'retinotopy-templates 11/22/2017 13:21:46';
[fwInfo] = getAnalysisFromFlywheel(theProject,analysisLabel,bensonDir, 'verbose', true, 'searchDir', projectDir);

%% Benson bookkeeping file 
% freesurfer file
cmd = ['echo -e "' 'This is a bookkeeping file.\nDo not delete unless you delete the whole analysis directory." > ' ...
       fullfile(sessionDir,'neuropythy',[fwInfo.analysis_id '.txt'])];
system(cmd);

%% Clean Up 
% remove the unused files from the fmriprep download
cmd = ['rm -rf ' fmriprepDir]; 
system(cmd)
