function [] = averageMaps(subjIds, mapOfInterest, varargin)
% Creates a scatter plot of a parameter of interest and eccentricity
%
% Syntax:
%   [figHndl] = scatterMapWithEcc(analysisParams,mapOfInterest, hemi)
%
% Description:
%    This function takes in a cell array of packets and returns one packet
%    with a the fields being a concatenation of all the input packets.
%
% Inputs:
%    subjIds                    - Cell array of subject IDs
%    mapOfInterest              - String with map name. Either:
%                                 'minorAxis'
%                                 'angle'
%                                 'amplitude'
%                                 'semi'
%                                 'exponent'
%                                 'rSqaured'
% Outputs:
%    figHndl                    - Figure Handle
%
% Optional key/value pairs:
%    saveFigs                   - Logical flag to save the figure
%    mapCoverage                - Nickname for the map vertices fit. This
%                                 is the folder name in the subject folder
%                                 in surfaceMaps. ('V1', 'EVC',
%                                 'wholebrain')

% MAB 03/10/20 Wrote it.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('subjIds',@iscell);
p.addRequired('mapOfInterest',@ischar);
p.addParameter('saveFigs',true,@islogical)
p.addParameter('dotColor',[1 0.5 0.5],@isvector)
p.addParameter('mapCoverage','EVC',@ischar)
p.parse(subjIds,mapOfInterest,varargin{:});


for ii = 1:length(subjIds)
    % load subject params
    analysisParams = getSubjectParams(subjIds{ii});
    
    %% LOAD THE PARAMTER MAP
    % Path the subject map data
    dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
    mapReadPath    = fullfile(dropBoxPath,'surfaceMaps',analysisParams.expSubjID,p.Results.mapCoverage);
    
    switch mapOfInterest
        case 'minorAxis'
            subjMapName = fullfile(mapReadPath,['minorAxisMap_', analysisParams.sessionNickname '_allAreas.dscalar.nii']);
        case 'angle'
            subjMapName = fullfile(mapReadPath,['angleMap_', analysisParams.sessionNickname '_allAreas.dscalar.nii']);
        case 'amplitude'
            subjMapName = fullfile(mapReadPath,['nlAmpMap_', analysisParams.sessionNickname '_allAreas.dscalar.nii']);
        case 'semi'
            subjMapName = fullfile(mapReadPath,['nlSemiMap_', analysisParams.sessionNickname '_allAreas.dscalar.nii']);
        case 'exponent'
            subjMapName = fullfile(mapReadPath,['nlExpMap_', analysisParams.sessionNickname '_allAreas.dscalar.nii']);
        case 'rSqaured'
            subjMapName = fullfile(mapReadPath,['rSquaredMapQcm_', analysisParams.sessionNickname '_allAreas.dscalar.nii']);
        case 'rSqauredDiff'
            subjMapName = fullfile(mapReadPath,['rSquaredDiffMap_', analysisParams.sessionNickname '.dscalar.nii']);
        case 'wholeBrainRsquared'
            subjMapName = fullfile(mapReadPath,['GLM_', analysisParams.sessionNickname '_wholeBrain.dscalar.nii']);
    end
    % Load the paramter map
    mapVals(:,ii) = loadCIFTI(subjMapName);
    
end

% Average the map of interest
if strcmp(mapOfInterest,'angle')
    angles = deg2rad(2*mapVals)';
    % Obtain the circular mean
    cMean = circ_mean(angles);
    % Convert back to deg
    averageMap = rad2deg(cMean)'./2;
else
    averageMap = mean(mapVals,2);
end

% Template File
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');
% Save path
mapSavePath =  fullfile(dropBoxPath,'surfaceMaps','averageSub');
switch mapOfInterest
    case 'minorAxis'
        outMapName = fullfile(mapSavePath,['minorAxisMap_Average_allAreas.dscalar.nii']);
    case 'angle'
        outMapName = fullfile(mapSavePath,['angleMap_Average_allAreas.dscalar.nii']);
    case 'amplitude'
        outMapName = fullfile(mapSavePath,['nlAmpMap_Average_allAreas.dscalar.nii']);
    case 'semi'
        outMapName = fullfile(mapSavePath,['nlSemiMap_Average_allAreas.dscalar.nii']);
    case 'exponent'
        outMapName = fullfile(mapSavePath,['nlExpMap_Average_allAreas.dscalar.nii']);
    case 'rSqaured'
        outMapName = fullfile(mapSavePath,['rSquaredMapQcm_Average_allAreas.dscalar.nii']);
    case 'rSqauredDiff'
        outMapName = fullfile(mapSavePath,['rSquaredDiffMap_Average_allAreas.dscalar.nii']);
    case 'wholeBrainRsquared'
        outMapName = fullfile(mapSavePath,['glm_r_squared_Average_wholebrain.dscalar.nii']);
end

makeWholeBrainMap(averageMap', [], templateFile, outMapName)

end
