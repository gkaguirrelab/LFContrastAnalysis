%% This is a demo to help us fit the QCM directly to the time course
%

%% Initialize
clear;

% Get subject specific params: 'LZ23', 'KAS25', 'AP26'
analysisParams = getSubjectParams('AP26');

% set the preprocessing method that was used to ananlyze the data.
analysisParams.preproc = 'hcp';

% turn on or off plotting
analysisParams.showPlots = true;

% Set the option to use simulated data from known parameters
analysisParams.analysisSimulate = false;
% Set which model to use to generate the
analysisParams.simulationMethod = 'QCM'; % 'QCM' or 'IAMP'

analysisParams.HRF = generateHRFKernel(6,12,10,analysisParams.timebase*1000);

%% Get the time series
% from simulation
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
    % From the data
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

%% Run the IAMP/QCM models
%
% Fit IAMP
%
% Fit IAMP to each constructed packet and create packetPocket cell array of
% all the fit packets.
%     packetPocket - Meta data of packePocket contains the direction/contrast form of the same packet.
%     iampOBJ - the tfe IAMP object
%     iampParams - cell array of iampParams for each object
%
% NOTE: Each session gets its own row in the packet pocket.  May want to sweep
% back at some point and match conventions in analysis params to this, for
% example by making the various cell arrays columns rather than rows to
% match.  Similarly with LMVectorAngles vector, which could turn into a
% matrix.
[analysisParams, iampTimeCoursePacketPocket, iampOBJ, iampParams, iampResponses, rawTC] = fit_IAMP(analysisParams,fullCleanData);

iampTimeCoursePacketPocket = {iampTimeCoursePacketPocket{1,:},iampTimeCoursePacketPocket{2,:}};

for ii = 1:length(iampTimeCoursePacketPocket)
    y(ii,:) = highpass(iampTimeCoursePacketPocket{ii}.response.values ,5/288,1/.8);
end
figure
counter =0;
for ii = 1:length(iampTimeCoursePacketPocket)
    if ii == 11
        figure
        counter = 0;
    end
    counter = counter +1;
    subplot(3,4,counter)
    hold on
    plot(iampTimeCoursePacketPocket{ii}.response.timebase,iampTimeCoursePacketPocket{ii}.response.values,'k')
    plot( iampTimeCoursePacketPocket{ii}.response.timebase, y(ii,:),'r')
    
    legend('original','filtered')
end


T = 0.8;
Fs = 1/T;
L = 288;
figure
counter =0;
for ii = 1:length(iampTimeCoursePacketPocket)
    if ii == 11
        figure
        counter = 0;
    end
    counter = counter +1;
    subplot(3,4,counter)
    hold on
    
    
    Y = fft(iampTimeCoursePacketPocket{ii}.response.values);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    Yf = fft(y(ii,:));
    P2f = abs(Yf/L);
    P1f = P2f(1:L/2+1);
    P1f(2:end-1) = 2*P1f(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    plot(f,P1,'k')
    plot(f,P1f,'r')
    xlabel('f (Hz)')
    ylabel('|fft|')
    legend('original','filtered')
    
end

