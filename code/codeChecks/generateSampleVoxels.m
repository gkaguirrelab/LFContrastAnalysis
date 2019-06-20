function [params theVoxels] = generateSampleVoxels(betaWeights,numDirections,numContrast,numVoxels, varargin)
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
%    noiseLevel                - The percent of white noise to be added, 0
%                                = no noise added (default = 0)
%    realExpParams             - BLock timing based on a real experimental
%                                run
 
% MAB 06/10/19


p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('betaWeights',@isvector);
p.addRequired('numDirections',@isnumeric);
p.addRequired('numContrast',@isnumeric);
p.addRequired('numPackets',@isnumeric);
p.addParameter('blockLength',15,@isnumeric);
p.addParameter('deltaT',800,@isnumeric);
p.addParameter('baselineCondNum',6,@isnumeric);
p.addParameter('noiseLevel',0,@isnumeric);
p.addParameter('realExpParams',[],@ismatrix);

p.parse(betaWeights,numDirections,numContrast,numVoxels,varargin{:});

blockLength = p.Results.blockLength;
deltaT      = p.Results.deltaT;
baselineCondNum = p.Results.baselineCondNum;
noiseLevel = p.Results.noiseLevel;

%% Step 1: Generate Stimuli
contrastCond = repmat(1:numContrast,[1,numDirections]);
directions = repelem(1:numDirections,numContrast);

randIndx = randperm(length(contrastCond));
dirContrast = [contrastCond;directions];
randDirContrast = dirContrast(:,randIndx);

start = 1:blockLength:length(contrastCond)*blockLength;
stop = blockLength:blockLength:length(contrastCond)*blockLength;

if isempty(p.Results.realExpParams)
    expParams = [start;stop;randDirContrast]';
else
    expParams = p.Results.realExpParams;
end

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

%% Step 3: Generate voxel responses
voxelResps = repmat(modelResponseStruct.values,[numVoxels,1]);

if noiseLevel ~=  0
    r = noiseLevel + ((-1.*noiseLevel)-noiseLevel).*rand(size(voxelResps));
    theVoxels = voxelResps + (voxelResps.*r);
else
    theVoxels = voxelResps;
end




