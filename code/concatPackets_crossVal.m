function [theTestPackets, theTrainPackets,leaveOutPairs] = concatPackets_crossVal(analysisParams, packetPocket, varargin)
% Takes in the packets and concatenates all the fields.
%
% Syntax:
%   [analysisParams, theFullPacket] = concatPackets(analysisParams, packetPocket)
%
% Description:
%    This function takes in a cell array of packets and returns one packet
%    with a the fields being a concatenation of all the input packets.
%
% Inputs:
%    analysisParams             - Struct of important information for the
%                                 analysis
%    packetPocket               - Cell array of packets (one packet per
%                                 run)
%
% Outputs:
%    analysisParams             - Returns analysisParams with any updates
%    theFullPacket              - The concatenated packet
%
% Optional key/value pairs:
%    bootstrap                  - Logical. If true, will randomly sample
%                                 with replacement the input packets.

% MAB 03/10/20 Wrote it.

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addRequired('packetPocket',@iscell);

p.parse(analysisParams,packetPocket,varargin{:});

%% Create the full timebase
trainTimePoints = analysisParams.expLengthTR*((analysisParams.numSessions*analysisParams.numAcquisitions)-2);
testTimePoints  = analysisParams.expLengthTR*((analysisParams.numSessions*analysisParams.numAcquisitions)-18);

trainTimeBase = linspace(0,(trainTimePoints-1)*analysisParams.TR,trainTimePoints)*1000;
testTimeBase  = linspace(0,(testTimePoints-1)*analysisParams.TR,testTimePoints)*1000;

%% Initialize the train packet
theTrainPacket.response.values   = [];
theTrainPacket.response.timebase = trainTimeBase;

% The stimulus
theTrainPacket.stimulus.values   = [];
theTrainPacket.stimulus.timebase = trainTimeBase;

% The kernel
kernelVec = zeros(size(trainTimeBase));
kernelVec(1:length(analysisParams.HRF.values)) = analysisParams.HRF.values;
theTrainPacket.kernel.values = kernelVec;
theTrainPacket.kernel.timebase = trainTimeBase;

% The metaData (this is the constrast and directions)
theTrainPacket.metaData.stimDirections = [];
theTrainPacket.metaData.stimContrasts  = [];
theTrainPacket.metaData.lmsContrast    = [];

%% Initialize the train packet
theTestPacket.response.values   = [];
theTestPacket.response.timebase = testTimeBase;

% The stimulus
theTestPacket.stimulus.values   = [];
theTestPacket.stimulus.timebase = testTimeBase;

% The kernel
kernelVec = zeros(size(testTimeBase));
kernelVec(1:length(analysisParams.HRF.values)) = analysisParams.HRF.values;
theTestPacket.kernel.values = kernelVec;
theTestPacket.kernel.timebase = testTimeBase;

% The metaData (this is the constrast and directions)
theTestPacket.metaData.stimDirections = [];
theTestPacket.metaData.stimContrasts  = [];
theTestPacket.metaData.lmsContrast    = [];

%% Put stuff in the cross val packets
lst = 11:20;
leaveOutPairs = [randperm(10)',lst(randperm(10))'];


for ii = 1:size(leaveOutPairs,1)
    
    trainOrder = 1:20;
    trainOrder([leaveOutPairs(ii,1),leaveOutPairs(ii,2)]) = [];
    
    testOrder = [leaveOutPairs(ii,1),leaveOutPairs(ii,2)];
    
    theTrainPackets{ii} = theTrainPacket;
    theTestPackets{ii}  = theTestPacket;
    
    %assemble train packet
    for jj = 1:length(trainOrder)
        % The Response
        theTrainPackets{ii}.response.values         = [theTrainPackets{ii}.response.values packetPocket{trainOrder(jj)}.response.values];
        % The Stimulus
        theTrainPackets{ii}.stimulus.values         = [theTrainPackets{ii}.stimulus.values packetPocket{trainOrder(jj)}.stimulus.values];
        % The metaData (this is the constrast and directions)
        theTrainPackets{ii}.metaData.stimDirections = [theTrainPackets{ii}.metaData.stimDirections packetPocket{trainOrder(jj)}.metaData.stimDirections];
        theTrainPackets{ii}.metaData.stimContrasts  = [theTrainPackets{ii}.metaData.stimContrasts  packetPocket{trainOrder(jj)}.metaData.stimContrasts];
        theTrainPackets{ii}.metaData.lmsContrast    = [theTrainPackets{ii}.metaData.lmsContrast packetPocket{trainOrder(jj)}.metaData.lmsContrast];
    end
    
    for kk = 1:length(testOrder)
        % The Response
        theTestPackets{ii}.response.values         = [theTestPackets{ii}.response.values packetPocket{testOrder(kk)}.response.values];
        % The Stimulus
        theTestPackets{ii}.stimulus.values         = [theTestPackets{ii}.stimulus.values packetPocket{testOrder(kk)}.stimulus.values];
        % The metaData (this is the constrast and directions)
        theTestPackets{ii}.metaData.stimDirections = [theTestPackets{ii}.metaData.stimDirections packetPocket{testOrder(kk)}.metaData.stimDirections];
        theTestPackets{ii}.metaData.stimContrasts  = [theTestPackets{ii}.metaData.stimContrasts  packetPocket{testOrder(kk)}.metaData.stimContrasts];
        theTestPackets{ii}.metaData.lmsContrast    = [theTestPackets{ii}.metaData.lmsContrast packetPocket{testOrder(kk)}.metaData.lmsContrast];
    end
    
    
end

end









