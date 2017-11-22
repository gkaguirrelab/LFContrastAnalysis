% ProcessDataFromFlywheelExample  Do some things with data downloaded from flywheel
%
% Description:
%   Assume we have brought down the output of fmriprep from flywheel, and
%   now do something good on it.
%
% See also:
%   GetDataFromFlywheelExample
%

% History
%  11/22/17  dhb, gka, mab  Tomorrow we give thanks, but today we code.

%% Clear
clear; close all;

%% Define our project
theProject = 'LFContrast';

%% Open flywheel object
fw = Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

%% Get some downloaded data
%
% A number of things hardcoded for now

% The label for the analysis we're going to look at. 
analysis_label = 'fmriprep 10/26/2017 22:17:09';
sesDirMagicName = 'ses-201709191435';

% Flywheel tree info, related to what gets unzipped
flywheelTree = 'flywheel/v0/output/';
fmriprepTree = 'fmriprep_output/fmriprep';

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

%% Find the output of fmri prep
%
% This code, in a somewhat hard coded manner, gets us to where the fmriprep
% output has been unzipped.  We cd to that directory, and in there are anat
% and func directories.
for ii = 1:numel(results)
    file_name = results(ii).file.name;
    output_name = fullfile(out_dir, file_name);
    
    session_id = results(ii).session.x_id;
    analysis_id = results(ii).analysis.x_id;
    
    % Figure out more about where things are in the tree.
    [~,body,ext] = fileparts(file_name);
    switch (ext)
        case '.html'
            % The name of the html file seems to be the name of the
            % directory that we unzip into.
            outputZipDir = body;
        case '.zip'
    end
end

% Assemble output directory and cd to there
fmriprepOutputDir = fullfile(out_dir,flywheelTree,fmriprepTree,outputZipDir,sesDirMagicName,[]);
curDir = pwd;
cd(fmriprepOutputDir);

%% Go get a functional files
%
% There should be some way to find out what these all are.  We can either
% do a dir command locally, or else perhaps there is a way to queyry the
% analysis up on flywheel.
niftiFunctionalFile = 'sub-HEROgka1_ses-201709191435_task-tfMRILFContrastPA_run-3_bold_space-MNI152NLin2009cAsym_preproc.nii.gz';

% Cd into the functional directory and load in a NIFTI file
cd('func');
fprintf('Loading in functional NIFTI file %s\n',niftiFunctionalFile);
funcData = load_nifti(niftiFunctionalFile);
fprintf('Read in a volume with %d elements\n',numel(funcData.vol));

%% Plot a time series
theVoxelRow = 45;
theVoxelCol = 17;
theVoxelSlice = 36;
timeSeries = squeeze(funcData.vol(theVoxelRow,theVoxelCol,theVoxelSlice,:));
figure; clf;
plot(timeSeries);

%% Plot slices
%
% timePoint = 8;
% for ss = 1:size(funcData.vol,3)
%     theSlice = funcData.vol(:,:,ss,timePoint);
%     figure(1); clf;
%     imagesc(theSlice);
%     drawnow;
% end







