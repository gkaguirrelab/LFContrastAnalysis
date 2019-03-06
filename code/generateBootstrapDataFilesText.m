function [] = generateBootstrapDataFilesText(analysisParams,sampleMatrix)
% This function writes a text file containing the random draws of the bootstap.
%
% Syntax:
%   [] = generateBootstrapDataFilesText(analysisParams,sampleMatrix);
%
% Description:
%    This function creates a text file containing the random draws of the
%    bootstrap in order for the analysis code to know the stim order for
%    each draw. This txt is named 'bootstrapDataFiles.txt' and is located
%    in the same folder as the dataFiles.txt for each session. This test
%    file is intended to be overwritten with each iteration of the bootstrap. 
%
% Inputs:
%    analysisParams            - Analysis parameter stuct set in analyzeLFContrast (Struct)
%    sampleMatrix              - A martix that is of size runs by sessions 
%                                that contains the random draw order of the
%                                bootstrap. Each collumn is a session that
%                                has the sampling with replacememnt if the
%                                bootstrap. 
%
% Outputs:
%    none
%
% Optional key/value pairs:
%    none

% MAB 03/04/19

for ii = 1:length(analysisParams.sessionNumber)
    
    %% Get file name
    bootstrapTrialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.sessionFolderName{ii},'experimentFiles','bootstrapDataFiles.txt');
    
    %% Open file
    fileID = fopen(bootstrapTrialOrderFile, 'w');
    
    
    %% Write to file
    % this should eventually be change to something of the form so it makes less assumptions 
    % fprintf(fileID,'CRF_%s_scan%1.f.mat\n',analysisParams.sessionNumber{ii},[3, 3, 4, 5,6])
    
    fprintf(fileID,'CRF_session1_scan%1.f.mat\n',sampleMatrix(:,ii))
    
    %% Close the file
    fclose(fileID);
end