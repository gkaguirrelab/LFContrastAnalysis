% IAMP Parameter Recovery Test

betaWeights = [repmat(5:-1:1,1,4), 0]';
numDirections = 4;
numContrast = 6;
numPackets = 4;
thePackets = generateSamplePackets(betaWeights,numDirections,numContrast,numPackets);

% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

for ii = 1:numPackets
    % Set the number of instances.
    clear defaultParamsInfo
    defaultParamsInfo.nInstances = size(thePackets{ii}.stimulus.values,1);
    
    [paramsFit{ii},fVal(ii),IAMPResponses{ii}] = iampOBJ.fitResponse(thePackets{ii},'defaultParamsInfo', ...
        defaultParamsInfo, 'searchMethod','linearRegression');
end

avgIampParams = iampOBJ.averageParams(paramsFit);


%% Concat and fit
concatPacket = iampOBJ.concatenatePackets(thePackets);
clear defaultParamsInfo
defaultParamsInfo.nInstances = size(concatPacket.stimulus.values,1);

[paramsFitConcat,fValConcat,IAMPResponsesConcat] = iampOBJ.fitResponse(concatPacket,'defaultParamsInfo', ...
    defaultParamsInfo, 'searchMethod','linearRegression');


%% Plot

figure; hold on
plot(betaWeights,'r')
plot(avgIampParams.paramMainMatrix,'b--')
