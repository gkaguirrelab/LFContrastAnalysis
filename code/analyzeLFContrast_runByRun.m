function [modelResponseStructIAMP, modelResponseStructQCM, thePacketIAMP, thePacketQCM]  = analyzeLFContrast_runByRun(subjId, varargin)
% Takes in a time series and chops it up into n runs.
%
% Syntax:
%    [thePackets] = chopUpTimeCourse(timeCoursePacket,numChops,varargin)
%
% Description:
%    Takes in a time series and chops it up nto n runs of equal length.
%
% Inputs:
%    timeCoursePacket           - The time course packet to be chopped up
%                                 A struct with "timebase" and "values"
%                                 subfeilds.
%    numChops                   - Number of cut to be made
%
% Outputs:
%    choppedTC                  - A cell array of the original packet
%                                 chopped into nunChops smaller packets
%                                 with a "values" and "timebase
% Optional key/value pairs:
%    - none for now

% MAB 12/22/19 created it

p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('subjId',@isstr);
p.addParameter('showPlot',@islogical);
p.parse(subjId,varargin{:});

%Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams(subjId);

%set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

%turn on or off plotting
analysisParams.showPlots = true;

%Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;
%Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

%set the HRF
load(fullfile(getpref('LFContrastAnalysis','melaAnalysisPath'),'LFContrastAnalysis','subjectHRFs',analysisParams.expSubjID,[analysisParams.expSubjID '_eventGain_results.mat']));
xBase = zeros(1,analysisParams.expLengthTR);
xBase(1:length(results.hrf')) = results.hrf';
analysisParams.HRF.values = xBase;
analysisParams.HRF.timebase =   analysisParams.timebase*1000;

hrfAUC = trapz(analysisParams.HRF.timebase,analysisParams.HRF.values);
analysisParams.HRF.values = analysisParams.HRF.values ./hrfAUC;

%analysisParams.HRF2 = generateHRFKernel(6,12,10,analysisParams.timebase*1000);

%Get the cleaned time series
if analysisParams.analysisSimulate
    analysisParams.numAcquisitions = 10;
    analysisParams.numSessions = 2;
    switch analysisParams.simulationMethod
        case 'IAMP'
            betaWeights = [repmat(1:-1/5:1/5,1,4), 0]';
            numDirections = 4;
            numContrast = 6;
            numVoxels = 400;
            [params,fullCleanData] = simulateDataFromExpParams(analysisParams,betaWeights,numDirections,numContrast,numVoxels, 'linDetrending', false);
        case 'QCM'
            angle = -45;
            minorAxisRatio = 0.19;
            fullCleanData = simulateDataFromEllipseParams(analysisParams,angle,minorAxisRatio,'numVoxels',850,...
                'crfOffset', 0,'noiseSD',2, 'noiseInverseFrequencyPower', .1);
    end
else
    switch analysisParams.preproc
        case 'fmriprep'
            [fullCleanData, analysisParams] = getTimeCourse(analysisParams);
        case 'hcp'
            [fullCleanData, analysisParams] = getTimeCourse_hcp(analysisParams);
        otherwise
            error('Preprocessing method unknown')
    end
end

% Run the IAMP/QCM models

% Fit IAMP

% Fit IAMP to each constructed packet and create packetPocket cell array of
% all the fit packets.
%    packetPocket - Meta data of packePocket contains the direction/contrast form of the same packet.
%    iampOBJ - the tfe IAMP object
%    iampParams - cell array of iampParams for each object
%
% NOTE: Each session gets its own row in the packet pocket.  May want to sweep
% back at some point and match conventions in analysis params to this, for
% example by making the various cell arrays columns rather than rows to
% match.  Similarly with LMVectorAngles vector, which could turn into a
% matrix.
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);



% Get directon/contrast form of time course and IAMP crf packet pockets.
% 
% This conversion is possible because the IAMP packet pocket has meta data
% that we put there to allow exactly this conversion.  That meta data
% encapsulates the key things we need to know about the stimulus obtained
% from the analysis parameters.
directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket);
directionTimeCoursePacketPocket = {directionTimeCoursePacketPocket{1,:},directionTimeCoursePacketPocket{2,:}};

% This puts together pairs of acquistions from the two sessions, so that
% we have one IAMP fit for each pair.  We do this because to fit the
% quadratic model, we need data for all of the color directions together.
% 
% NOTE: This bit is very specific to the design of the experiment we are
% currently analyzing, and has to do specifically with the way color
% directions were studied across acquisitions and sessions.
% 
% 
% ###### FIX ###################
% remove subraction of the baseline
% ##############################
for ii = 1:analysisParams.numAcquisitions
    [concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','averageBaseline');
    %[concatParams{ii},concatBaselineShift(:,ii)] = iampOBJ.concatenateParams(iampParams(:,ii),'baselineMethod','makeBaselineZero');
end

medianIampParams = iampOBJ.medianParams(concatParams);

directionCrfMeanPacket = makeDirectionCrfPacketPocket(analysisParams,medianIampParams);

% Fit the CRF with the QCM -- { } is because this expects a cell
[qcmCrfMeanOBJ,qcmParams] = fitDirectionModel(analysisParams, 'qcmFit', directionTimeCoursePacketPocket,'fitErrorScalar',1000);
qcmMedianParams = qcmCrfMeanOBJ.medianParams(qcmParams)
qcmCrfMeanOBJ.paramPrint(qcmMedianParams)
end
