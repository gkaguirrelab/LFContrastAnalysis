function voxelTimeSeries = simulateDataFromEllipseParams(analysisParams,angle,minorAxisRatio, varargin)
% Generates simulated voxel timecourses based on betaweights and
%
% Syntax:
%    thePackets = generateSamplePackets(betaWeights,numDirections,numContrast,numPackets)
%
% Description:
%    This function takes in a vector of beta weights (one for each
%    regressor) and returns a random time course repsonse packet.
%
% Inputs:
%    analysisParams            - Beta weight values for each condtion
%    betaWeights               -  Number of directions
%    numDirections             - Number of contrast conditions
%    numContrast               - Number of packets to be generated
%    numVoxels                 -
% Outputs:
%    thePackets                - Packets of simulated time courses
% Optional key/value pairs:
%    baseline                  - Baseline added to each run
%    linDetrending             -

% MAB 06/10/19

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isvector);
p.addRequired('angle',@isnumeric);
p.addRequired('minorAxisRatio',@isnumeric);
p.addParameter('crfAmp',8,@isnumeric);
p.addParameter('crfExponent',1.2,@isnumeric);
p.addParameter('crfSemi',0.9,@isnumeric);
p.addParameter('expFalloff',5,@isnumeric);
p.addParameter('crfOffset',-.02,@isnumeric);
p.addParameter('numVoxels',900,@isnumeric);
p.addParameter('addNoise',true,@islogical);
p.addParameter('noiseLevel',0.2,@isnumeric);

p.parse(analysisParams, angle, minorAxisRatio, varargin{:});

% Initialize the fit object for the QCM
fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);

%% make the params struct
params.Qvec        = [minorAxisRatio, angle];
params.crfAmp      = p.Results.crfAmp;
params.crfExponent = p.Results.crfExponent;
params.crfSemi     = p.Results.crfSemi;
params.expFalloff  = p.Results.expFalloff;
params.crfOffset   = p.Results.crfOffset;

for sessionNum = 1:analysisParams.numSessions
    trialOrderDir  = fullfile(getpref(analysisParams.projectName,'projectPath'), analysisParams.projectNickname, 'DataFiles', analysisParams.expSubjID,analysisParams.sessionDate{sessionNum},analysisParams.sessionNumber{sessionNum});
    trialOrderFile = fullfile(getpref(analysisParams.projectName,'melaAnalysisPath'),'LFContrastAnalysis',analysisParams.sessionFolderName{sessionNum},'experimentFiles','dataFiles.txt');
    trialOrderFiles = textFile2cell(trialOrderFile);
    
    % get then number of direction over all sessiosn
    analysisParams.numDirections = size(unique(analysisParams.directionCoding','rows')',2);
    
    % Get the directions and contrasts for corresponding session
    if analysisParams.numDirPerSession == analysisParams.numDirections
        directionCoding  = analysisParams.directionCoding;
        maxContrast = analysisParams.maxContrastPerDir;
    elseif analysisParams.numDirPerSession < size(unique(analysisParams.directionCoding','rows')',2)
        sPos = 1+ analysisParams.numDirPerSession*(sessionNum-1);
        ePos = (1+ analysisParams.numDirPerSession*(sessionNum-1)) + (analysisParams.numDirPerSession-1);
        directionCoding = analysisParams.directionCoding(:,sPos:ePos);
        maxContrast = analysisParams.maxContrastPerDir(sPos:ePos);
    end
    
    
    for jj = 1:analysisParams.numAcquisitions
        
        %% Create a packet with the stimulus.values field populated by the stimulus design matrix (conditions x timepoints)
        % load relevent experiment info
        dataParamFile = fullfile(trialOrderDir,trialOrderFiles{jj});
        load(dataParamFile);
        
        % Get a matrix of block onset and offset TRs and contiditons
        expParams = getExpParams(dataParamFile,analysisParams.TR,'hrfOffset', false, 'stripInitialTRs', false);
        
        % make timebase
        totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
        deltaT = analysisParams.TR*1000;
        timebase = linspace(0,totalTime-deltaT,totalTime/deltaT);
        
        theKernel = generateHRFKernel(6,12,10,timebase);
        
        % make stimulus values for QCM
        contrastCoding = [analysisParams.contrastCoding, 0];
        LMSContrastMat = LMSContrastValuesFromParams(expParams,contrastCoding,directionCoding,maxContrast,totalTime,deltaT);
        directionPrecision = 4;
        indDirectionDirections = round(directionCoding(1:analysisParams.theDimension,:),directionPrecision);
        LMSContrastMat(3,:) = [];
        [stimDirections,stimContrasts] = tfeQCMStimuliToDirectionsContrasts(LMSContrastMat, ...
            'zeroContrastDirection',indDirectionDirections(:,1),'precision',directionPrecision);
        
        stimulusStruct.values   = [stimDirections;stimContrasts];
        stimulusStruct.timebase = timebase;
        
        % Generate time course prediction
        modelPreds = fitOBJ.computeResponse(params,stimulusStruct,theKernel);
        
        
        voxelTimeSeries(:,:,jj + analysisParams.numAcquisitions*(sessionNum-1)) = repmat(modelPreds.values,[p.Results.numVoxels,1]);
        
        
    end
    if p.Results.addNoise
        maxSignal = max(voxelTimeSeries(:));
        minSignal = min(voxelTimeSeries(:));
        theNoise  = p.Results.noiseLevel .* (minSignal + (maxSignal-minSignal).*rand(size(voxelTimeSeries)));
        voxelTimeSeries = voxelTimeSeries + theNoise;
    end
end








