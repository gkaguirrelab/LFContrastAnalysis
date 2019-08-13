function [movenentRegressors] = hcpReadMovementFiles(txtFile,rws,cols,varargin)
% Takes in the analysis params struct and returns a voxel by timepoint by
% aquisistion matrix for all the runs found fro the specicied session(s)
%
% Syntax:
%   [fullCleanData, analysisParams, voxelIndex] = getTimeCourse(analysisParams)
%
% Description:
%    This function takes in a struct that is specified in analyzeLFContrast.m
%    and returns a voxel by timepoint by aquisition matrix for all the
%    aquisitions specified in the text files housed in mela_analysis that
%    describe the session(s). f mutliple sessions, they will be
%    concatenated in the 3rd dimension
%
% Inputs:
%    analysisParams    - Stuct contianing relevenat info to the session that
%                        is defined in analyzeLFContrast.m. (string)
%
% Outputs:
%    fullCleanData       - The voxel by timepoint by aquisition matrix
%    analysisParams      - the input analysis params updated with the
%                          number of aquistidiond found per session
%    voxelIndex          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    selectCols
%    plotMovement

% MAB 09/09/18

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('txtFile',@isstr);
p.addRequired('rws',@isnumeric);
p.addRequired('cols',@isnumeric);
p.addParameter('selectCols',[],@isnumeric);
p.addParameter('plotMovement',false,@islogical);

p.parse(txtFile,rws,cols,varargin{:})

fileID = fopen(txtFile,'r');
formatSpec = '%f';
textVector = fscanf(fileID,formatSpec);
fclose(fileID);

movementRegressors = reshape(textVector,[cols,rws])';

if ~isempty(p.Results.selectCols)
    movenentRegressors = movementRegressors(:,p.Results.selectCols);
end

if p.Results.plotMovement
    
    colLabels = {'trans x (mm)', 'trans y (mm)', 'trans z (mm)', 'rot x (deg)' ,'rot y (deg)', 'rot z (deg)'};
    
    figure;
    plot(movementRegressors);
    
    
    % put info
    ylabel('Movement (mm or deg)')
    xlabel('Frames')
    
    
    if isempty(p.Results.selectCols)
        labelIndx = 1:lenth(colLabels);
    else
        labelIndx = p.Results.selectCols;
    end
    
    legend(colLabels(labelIndx))
end

end

