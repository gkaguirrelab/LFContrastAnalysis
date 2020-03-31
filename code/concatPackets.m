function [analysisParams, theFullPacket] = concatPackets(analysisParams, packetPocket, varargin)
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
p.addParameter('bootstrap',false,@islogical)

p.parse(analysisParams,packetPocket,varargin{:});

% Creat the full timebase
numTimePoints = analysisParams.expLengthTR*analysisParams.numSessions*analysisParams.numAcquisitions;
timebase = linspace(0,(numTimePoints-1)*analysisParams.TR,numTimePoints)*1000;

% Initialize the packet
theFullPacket.response.values   = [];
theFullPacket.response.timebase = timebase;

% The stimulus
theFullPacket.stimulus.values   = [];
theFullPacket.stimulus.timebase = timebase;

% The kernel
kernelVec = zeros(size(timebase));
kernelVec(1:length(analysisParams.HRF.values)) = analysisParams.HRF.values;
theFullPacket.kernel.values = kernelVec;
theFullPacket.kernel.timebase = timebase;

% The metaData (this is the constrast and directions)
theFullPacket.metaData.stimDirections = [];
theFullPacket.metaData.stimContrasts  = [];
theFullPacket.metaData.lmsContrast    = [];

 if p.Results.bootstrap == true
        runOrder = [ randi([1 10],1,10) ,  randi([11 20],1,10)];
    else 
        runOrder = 1:20;
 end
    
for ii = 1:length(packetPocket)

    % The Response
    theFullPacket.response.values   = [theFullPacket.response.values packetPocket{runOrder(ii)}.response.values];
    % The Stimulus
    theFullPacket.stimulus.values   = [theFullPacket.stimulus.values packetPocket{runOrder(ii)}.stimulus.values];
    % The metaData (this is the constrast and directions)
    theFullPacket.metaData.stimDirections = [theFullPacket.metaData.stimDirections packetPocket{runOrder(ii)}.metaData.stimDirections];
    theFullPacket.metaData.stimContrasts  = [theFullPacket.metaData.stimContrasts  packetPocket{runOrder(ii)}.metaData.stimContrasts];
    theFullPacket.metaData.lmsContrast    = [theFullPacket.metaData.lmsContrast packetPocket{runOrder(ii)}.metaData.lmsContrast];
end


end




