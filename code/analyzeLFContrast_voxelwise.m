function [] = analyzeLFContrast_voxelwise(subjId)

display(['STARTING - Making Maps: ',subjId])
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);
sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
mapSavePath    = fullfile(dropBoxPath,'surfaceMaps',analysisParams.expSubjID);
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');

%% Set up stuff
% Initialize the maps
qcmParamMap = zeros(91282,7);
meanRSquaredQCMMap = zeros(91282,1);
stdRSquaredQCMMap = zeros(91282,1);
meanRSquaredIAMPMap = zeros(91282,1);
stdRSquaredIAMPMap = zeros(91282,1);

% Create the fit object 
fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);

%% Get the cleaned time series
[fullCleanData, analysisParams, voxelIndex] = getTimeCourse_hcp(analysisParams);
% reshape the data to voxels x time point(all 20 runs)
timeCourses = [];
for ii = 1:size(fullCleanData,3)
    timeCourses = [timeCourses,fullCleanData(:,:,ii)];
end

%% load the HRF
[analysisParams] = loadHRF(analysisParams);


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
    
    % Fit the time course with the QCM -- { } is because this expects a cell
    [qcmTcOBJ,qcmParams] = fitDirectionModel(analysisParams, 'qcmFit', {theFullPacket},'fitErrorScalar',1000,'talkToMe',false);
    
    % Compute the response 
    theModelPred = fitOBJ.computeResponse(qcmParams{1},theFullPacket.stimulus,theFullPacket.kernel);
    qcmCorrVec =  [theFullPacket.response.values',theModelPred.values'];
    qcmCorrVals = corrcoef(qcmCorrVec(:,1),qcmCorrVec(:,2),'rows','complete').^2;

    % Put values back in map
    qcmParamMap(voxelIndex(ii),:)         = fitOBJ.paramsToVec(qcmParams{1});
    rSquaredQcmMap(voxelIndex(ii),:)      = qcmCorrVals(1,2);
end


% write out minor axis ratio
ciftiVec = qcmParamMap(:,1);
mapName        = fullfile(mapSavePath,['minorAxisMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = qcmParamMap(:,2);
mapName        = fullfile(mapSavePath,['angleMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = qcmParamMap(:,3);
mapName        = fullfile(mapSavePath,['nlAmpMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = qcmParamMap(:,4);
mapName        = fullfile(mapSavePath,['nlSemiMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out minor axis ratio
ciftiVec = qcmParamMap(:,5);
mapName        = fullfile(mapSavePath,['nlExpMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out mean R squared map
mapName        = fullfile(mapSavePath,['rSquaredMapQcm_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(rSquaredQcmMap', [], templateFile, mapName)
display(['COMPLETED: ',subjId])
end