function [] = analyzeLFContrast_voxelwise(subjId,varargin)
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('subjID',@ischar);
p.addParameter('roiType','V1',@ischar);
p.parse(subjId,varargin{:});


display(['STARTING - Making Voxelwise Maps: ',subjId])
% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId, 'roiType',p.Results.roiType);

sessionDir     = fullfile(getpref(analysisParams.projectName,'projectRootDir'),analysisParams.expSubjID);
dropBoxPath     = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),analysisParams.projectName);
mapSavePath    = fullfile(dropBoxPath,'surfaceMaps',analysisParams.expSubjID);
templateFile   = fullfile(dropBoxPath,'surfaceMaps','templates','template.dscalar.nii');

%% Set up stuff
% Initialize the maps
qcmParamMap = zeros(91282,7);
rSquaredQcmMap = zeros(91282,1);
rSquaredIampMap = zeros(91282,1);
rSquaredDiffMap = zeros(91282,1);
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

% the iamp
iampOBJ = tfeIAMP('verbosity','none');
defaultParamsInfo.nInstances = size(theConcatPacket.stimulus.values,1);

%% loop over voxels
for ii = 1:size(fullCleanData,1)
    
    % add the voxel time course to the packet
    theFullPacket.response.values = timeCourses(ii,:);
    theConcatPacket.response.values = timeCourses(ii,:);
    
    % Fit the time course with the QCM -- { } is because this expects a cell
    [qcmTcOBJ,qcmParams] = fitDirectionModel(analysisParams, 'qcmFit', {theFullPacket},'fitErrorScalar',1000,'talkToMe',false);
    
    [iampParams,fVal,iampResponses] = iampOBJ.fitResponse(theConcatPacket,...
        'defaultParamsInfo', defaultParamsInfo, 'searchMethod','linearRegression');
    % Compute the response
    theModelPred = qcmTcOBJ.computeResponse(qcmParams{1},theFullPacket.stimulus,theFullPacket.kernel);
    
    %% nan mask
    if any(isnan(theConcatPacket.response.values))
        nanMask =  find(isnan(theConcatPacket.response.values));
        theConcatPacket.response.values(nanMask) = [];
        theFullPacket.response.values(nanMask) =[];
        theModelPred.values(nanMask) = [];
        iampResponses.values(nanMask) = [];
        
        qcmCorrVals = corrcoef(theFullPacket.response.values',theModelPred.values').^2;
        
        iampCorrVals = corrcoef(theConcatPacket.response.values',iampResponses.values').^2;
        
    else
        
        qcmCorrVals = corrcoef(theFullPacket.response.values',theModelPred.values').^2;
        
        iampCorrVals = corrcoef(theConcatPacket.response.values',iampResponses.values').^2;
    end
    % Put values back in map
    qcmParamMap(voxelIndex(ii),:)         = qcmTcOBJ.paramsToVec(qcmParams{1});
    rSquaredQcmMap(voxelIndex(ii),:)      = qcmCorrVals(1,2);
    
    % Put values back in map
    rSquaredIampMap(voxelIndex(ii),:)      = iampCorrVals(1,2);
    
    rSquaredDiffMap(voxelIndex(ii),:)      = iampCorrVals(1,2)- qcmCorrVals(1,2);
end


% write out minor axis ratio
ciftiVec = qcmParamMap(:,1);
mapName        = fullfile(mapSavePath,['minorAxisMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out angle
ciftiVec = qcmParamMap(:,2);
mapName        = fullfile(mapSavePath,['angleMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out amp
ciftiVec = qcmParamMap(:,3);
mapName        = fullfile(mapSavePath,['nlAmpMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out semi
ciftiVec = qcmParamMap(:,4);
mapName        = fullfile(mapSavePath,['nlSemiMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out exp
ciftiVec = qcmParamMap(:,5);
mapName        = fullfile(mapSavePath,['nlExpMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out offset
ciftiVec = qcmParamMap(:,7);
mapName        = fullfile(mapSavePath,['nlOffset_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(ciftiVec', [], templateFile, mapName)

% write out mean QCM R squared map
mapName        = fullfile(mapSavePath,['rSquaredMapQcm_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(rSquaredQcmMap', [], templateFile, mapName)

% write out mean IAMP R squared map
mapName        = fullfile(mapSavePath,['rSquaredMapIamp_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(rSquaredIampMap', [], templateFile, mapName)

% write out mean R squared diff map
mapName        = fullfile(mapSavePath,['rSquaredDiffMap_', analysisParams.sessionNickname '.dscalar.nii']);
makeWholeBrainMap(rSquaredDiffMap', [], templateFile, mapName)

display(['COMPLETED: ',subjId])
end