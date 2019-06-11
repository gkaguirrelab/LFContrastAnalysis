%IAMP Parameter Recovery Test
% Demo script to test the IAMP parameter recovery

%% Set up parameters
%
% set the beta weights
betaWeights = [repmat(5:-1:1,1,4), 0]';

% number of directions
numDirections = 4;

% number of contrasts
numContrast = 6;

% number of pakets to generate
numPackets = 4;

%% Generate the packet
thePackets = generateSamplePackets(betaWeights,numDirections,numContrast,numPackets);

%% Fit the IAMP to seperate runs and average the beta weights
%
% Construct the model object
iampOBJ = tfeIAMP('verbosity','none');

% Loop over packets
for ii = 1:numPackets
    
    % Set the number of instances.
    clear defaultParamsInfo
    defaultParamsInfo.nInstances = size(thePackets{ii}.stimulus.values,1);
    
    % Use the IAMP Obj to fit the packet
    [paramsFit{ii},fVal(ii),IAMPResponses{ii}] = iampOBJ.fitResponse(thePackets{ii},'defaultParamsInfo', ...
        defaultParamsInfo, 'searchMethod','linearRegression');
end

avgIampParams = iampOBJ.averageParams(paramsFit);

%% Concatenate the packets and fit the IAMP 
%
% Concatenate the generated packets
concatPacket = iampOBJ.joinPackets(thePackets);

% Set the number of instances.
clear defaultParamsInfo
defaultParamsInfo.nInstances = size(concatPacket.stimulus.values,1);

% Use the IAMP Obj to fit the concatenated packet
[paramsFitConcat,fValConcat,IAMPResponsesConcat] = iampOBJ.fitResponse(concatPacket,'defaultParamsInfo', ...
    defaultParamsInfo, 'searchMethod','linearRegression');

%% Plot the recovered beta weights
figure; hold on
plot(betaWeights,'r')
plot(avgIampParams.paramMainMatrix,'b--')
plot(paramsFitConcat.paramMainMatrix,'g*')
legend('Original', 'Average', 'Concatenated')
