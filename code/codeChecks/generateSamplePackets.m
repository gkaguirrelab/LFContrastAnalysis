function thePackets = generateSamplePackets(betaWeights,numDirections,numContrast,numPackets, varargin)
% Generates simulated packets for the IAMP fitting routine for LFContrast
%
% Syntax:
%    thePackets = generateSamplePackets(betaWeights,numDirections,numContrast,numPackets)
%
% Description:
%    This function takes in a vector of beta weights (one for each
%    regressor) and returns a random time course repsonse packet. 
%
% Inputs:
%    betaWeights               - Beta weight values for each condtion 
%    numDirections             - Number of directions
%    numContrast               - Number of contrast conditions
%    numPackets                - Number of packets to be generated
% Outputs:
%    thePackets                - Packets of simulated time courses
% Optional key/value pairs:
%    blockLength               - The length of each block in TRs
%    deltaT                    - TR length in mS
%    baselineCondNum           - Contrast condition number corresponding to
%                                the baseline blocks

% MAB 06/10/19


p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('betaWeights',@isvector);
p.addRequired('numDirections',@isnumeric);
p.addRequired('numContrast',@isnumeric);
p.addRequired('numPackets',@isnumeric);
p.addParameter('blockLength',15,@isnumeric);
p.addParameter('deltaT',800,@isnumeric);
p.addParameter('baselineCondNum',6,@isnumeric);

p.parse(betaWeights,numDirections,numContrast,numPackets,varargin{:});

blockLength = p.Results.blockLength;
deltaT      = p.Results.deltaT;
baselineCondNum = p.Results.baselineCondNum;

for ii = 1:numPackets
    %% Step 1: Generate Stimuli
    contrastCond = repmat(1:numContrast,[1,numDirections]);
    directions = repelem(1:numDirections,numContrast);
    
    randIndx = randperm(length(contrastCond));
    dirContrast = [contrastCond;directions];
    randDirContrast = dirContrast(:,randIndx);
    
    start = 1:blockLength:length(contrastCond)*blockLength;
    stop = blockLength:blockLength:length(contrastCond)*blockLength;
    
    expParams = [start;stop;randDirContrast]';
    
    stimRegressor =  createRegressors(expParams,baselineCondNum,deltaT*stop(end),deltaT);
    
    params.paramNameCell = {'amplitude'};
    params.paramMainMatrix = betaWeights;
    params.matrixRows = size(betaWeights,1);
    params.matrixCols  = size(betaWeights,2);
    params.noiseSd = 0;
    
    
    %% Step 2: Generate Response
    stimulusStruct.values = stimRegressor;
    stimulusStruct.timebase = deltaT:deltaT:deltaT*stop(end);
    kernelStruct.timebase = deltaT:deltaT:5000;
    kernelStruct = generateHRFKernel(6,12,10,kernelStruct.timebase);
    
    % Construct the model object
    iampOBJ = tfeIAMP('verbosity','none');
    
    modelResponseStruct = iampOBJ.computeResponse(params,stimulusStruct,kernelStruct);
    
    %% Step 3: Genmerate Packet
    thePackets{ii}.response.values   = modelResponseStruct.values;
    thePackets{ii}.response.timebase = modelResponseStruct.timebase;
    % the stimulus
    thePackets{ii}.stimulus.timebase = stimulusStruct.timebase;
    thePackets{ii}.stimulus.values   = stimulusStruct.values;
    thePackets{ii}.stimulus.metaData = [];
    % the kernel
    thePackets{ii}.kernel = kernelStruct;
    % the meta data (this is the constrast and directions)
    thePackets{ii}.metaData.stimDirections = [];
    thePackets{ii}.metaData.stimContrasts  = [];
    thePackets{ii}.metaData.lmsContrast    = [];
    
end