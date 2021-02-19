function [] = analyzeLFContrast_wholeBrainGLM(subjId)

display(['STARTING - Making Maps: ',subjId])
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);
sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
mapSavePath    = fullfile(dropBoxPath,'surfaceMaps',analysisParams.expSubjID);
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');

%% load the HRF
[analysisParams] = loadHRF(analysisParams);

%% Set up stuff
% Initialize the maps
rSquaredGlmMap = zeros(91282,1);

% Create the fit object
iampOBJ = tfeIAMP('verbosity','none');

%% Get the cleaned time series
w = warning ('off','all'); % turn off waring for nuisance regression
[fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams,'wholeBrain',true);
warning(w)
% reshape the data to voxels x time point(all 20 runs)
timeCourses = [];
for ii = 1:size(fullCleanData,3)
    timeCourses = [timeCourses,fullCleanData(:,:,ii)];
end

%% Initialize the packet
% Get a packet for each run (1-20)
[analysisParams, iampTimeCoursePacketPocket] = generateRunPackets(analysisParams, fullCleanData);

% Generate a cell array of concat train packets and a corresponding cell
%  array of concat test packets
[analysisParams, theConcatPacket] = concatPackets(analysisParams, iampTimeCoursePacketPocket);
tmpPkt = makeDirectionTimeCoursePacketPocket({theConcatPacket});
theFullPacket = tmpPkt{1};
theFullPacket.response.values = [];

%% loop over voxels
for ii = 1:size(fullCleanData,1)
    
    % add the voxel time course to the packet
    theFullPacket.response.values = timeCourses(ii,:);
    
    % fit the IAMP model
    defaultParamsInfo.nInstances = size(theFullPacket.stimulus.values,1);
    [iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theFullPacket,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    
    % Compute the response
    
    glmCorrMat =  [theFullPacket.response.values',iampResponses.values'];
    glmCorrVals = corrcoef(glmCorrMat(:,1),glmCorrMat(:,2),'rows','complete').^2;
    
    % Put values back in map
    rSquaredGlmMap(voxelIndex(ii),:)      = glmCorrVals(1,2);
end


% write out GLM R^2 map
mapName        = fullfile(mapSavePath,['GLM_', analysisParams.sessionNickname '_wholeBrain.dscalar.nii']);
makeWholeBrainMap(rSquaredGlmMap', [], templateFile, mapName)

display(['COMPLETED: ',subjId])
end