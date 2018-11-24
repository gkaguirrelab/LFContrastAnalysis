% Test the QCM fitting and related matters
%
% Description:
%   This script synthesizes data for the QCM model and then fits it,
%   to ensure that we can get back what we put in.  Etc.
% 
%   Note that this is developed for the two-dimensional (ellipse) version.
% 
%   NOTE: When we get to ellipses, we need to remember to put a constraint into the
%   fitting that keeps the length of the third axis smaller than the second.
%

% History:
%   11/20/18  dhb, mab  Tuned this up; things are making sense.
%   1//24/18  dhb       A bit more tuning.  Fix row/col convention.

%% Initialize
clear; close all

%% Set up params
%
% Naka Ruston Params
NOOFFSET = false;
theDimension = 2;
Rmax   = 0.9;
sigma  = 0.7;
n      = 2.1;
if (NOOFFSET)
    offset = 0;
else
    offset = -0.1;
end

% Ellipse parameters
minorAxis = 0.3;
rotdeg    = 45;
ellParams = [1 minorAxis rotdeg];

%% Generate Stimulus:
%
% Random stimuli bounded between 1 and -1
RANDOM_STIMULI = false;
if (RANDOM_STIMULI)
    numStim = 300;
    stim = (2*rand(2,numStim) -1);
    
% Stimuli that vary along specified directions
else
    maxContrast = 0.6;
    contrastCoding = [0.0625, 0.125, 0.25, 0.5, 1];
    nContrastsPerDirection = length(contrastCoding);
    theDirections = [ [1 0]',  [1 1]',  [0 1]', [1 -1]' ];
    for ii = 1:size(theDirections,2)
        theDirections(:,ii) = theDirections(:,ii)/norm(theDirections(:,ii));
    end
    directionCoding = maxContrast*theDirections; 
    stim = kron(directionCoding',contrastCoding')';
    numStim = size(stim,2);
end

%% Generate response
% This does by hand what our compute response routine does, and is here
% just to expose that calculation and make sure it does what we think it
% should.
%
% Set up scale matrix
S = diag(ellParams(1:2));

% Set up rotation matrix
%
% Whether or not we tag on a transpose here defines 
% the convention for angle. Currently matched to the
% tfeQCM implementation and seems like the right convention
% for the two-dimensional case at least.
V = deg2rotm(ellParams(3))';  

% Get the Q matrix that takes stimulus and get a radius. 
A = S*V';
Q = A'*A;

% Get the radius
radius = diag(sqrt(stim'*Q*stim));

% Get the neural response values
R  = nakaRushton(radius,sigma,n,Rmax,offset);

%% Let's check that the QCM forward model gives the same responses.
%
% Construct the model object 
tfeResponseCheck = tfeQCM('verbosity','none','dimension',theDimension);
stimulusStruct.values = stim;
stimulusStruct.timebase = 1:numStim;
nTimeSamples = size(stimulusStruct.timebase,2);

% Set parameters and simulate responses
params1 = tfeResponseCheck.defaultParams;
params1.Qvec = [minorAxis rotdeg];
params1.crfAmp = Rmax;
params1.crfSemi = sigma;
params1.crfExponent = n;
params1.noiseSd = 0.01;
params1.offset = offset;
modelResponseStruct = tfeResponseCheck.computeResponse(params1,stimulusStruct,[],'AddNoise',false);
if (max(abs(R'-modelResponseStruct.values)) > 1e-15)
    error('Hand computation of QCM model does not match tfeQCM forward model');
end


%%  Use the tfeQCM to fit the stim/resp:
%
% Get the tfeQCM object
temporalFitQCM = tfeQCM('verbosity','none','dimension',theDimension);

% Set up the packet with the stimulus
stimulusStruct.values   = stim;
stimulusStruct.timebase = 1:length(stimulusStruct.values);

% Set up the packet with the response
thePacket.response.values = R';
thePacket.response.timebase = 1:length(thePacket.response.values);

% Construct a packet for the QCM to fit.
thePacket.stimulus = stimulusStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Fit the packet
if (NOOFFSET)
    defaultParamsInfo.noOffset = true;
else
    defaultParamsInfo.noOffset = false;
end
[paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket,'defaultParamsInfo',defaultParamsInfo);
fprintf('Model parameter from fits:\n');
temporalFitQCM.paramPrint(paramsQCMFit)

%%  Check that the fit recovers the responses we put in
% This is done by hand to match hand-coded method above,
% but could be done by a call through the tfeQCM routine.
qcmParams =[1 paramsQCMFit.Qvec];
S_qcm = diag(qcmParams(1:2));

% Set up rotation matrix
V_qcm = deg2rotm(qcmParams(3))';  

% Get the Q matrix that takes stim and get a radius. 
A_qcm = S_qcm*V_qcm';
Q_qcm = A_qcm'*A_qcm;

% Get the raduis
radiusQcm =  diag(sqrt(stim'*Q_qcm*stim));

% Get the predicted response values.  These should match reasonably well
% the responses we simulated above.  But note that there are some
% ambiguities in the parameterization of the QCM model, because of a +/- 90
% degree ambiguit about which is the major and which is the minor axis of
% the ellipse.
Rqcm  = nakaRushton(radiusQcm,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent,paramsQCMFit.crfAmp, paramsQCMFit.offset);
if (max(abs(R-Rqcm)/max(abs(R))) > 1e-2)
    error('Hand computation of QCM model does not match tfeQCM forward model');
end

%%  Plot simulated and predicted responses
figure; hold on
p1 = scatter3(stim(1,:),stim(2,:),radius,'b','*');
p2 = scatter3(stim(1,:),stim(2,:),radiusQcm,'r');
xlabel('L Contrast')
ylabel('M Contrast')
zlabel('Radius')
title('The radius as found by s''*Q*s = r');
legend([p1 p2], 'original', 'QCM recovered')

figure; hold on 
p3 = scatter3(stim(1,:),stim(2,:),R,'b','*');
p4 = scatter3(stim(1,:),stim(2,:),Rqcm,'r');
legend([p3 p4], 'original', 'QCM recovered')
xlabel('L Contrast')
ylabel('M Contrast')
zlabel('Response')
title('The response as found by passing the radius throught the naka rushton');
legend([p1 p2], 'original', 'QCM recovered')

%%  Check that Naka-Rushton funciton inverts
thresh = 0.3;
eqContrast = InvertNakaRushton([paramsQCMFit.crfAmp,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent],thresh);
circlePoints = eqContrast*UnitCircleGenerate(numStim);
[~,Ainv,Q] = EllipsoidMatricesGenerate([1 paramsQCMFit.Qvec],'dimension',2);
ellipsePoints = Ainv*circlePoints;
checkThresh = ComputeNakaRushton([paramsQCMFit.crfAmp,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent],diag(sqrt(ellipsePoints'*Q*ellipsePoints)));
if (any(abs(checkThresh-thresh) > 1e-10))
    error('Did not invert QCM model correctly');
end

%%  Find contrast in given direction that produces desired response
if (~RANDOM_STIMULI)
    maxResponseFactor = 3;
    whichDirection = 2;
    
    % Pull out this direction and simulated responses
    theDirection = theDirections(:,whichDirection);
    directionResponses = R(nContrastsPerDirection*(whichDirection-1)+1:nContrastsPerDirection*whichDirection);
    maxResponse = max(directionResponses);
    
    % Invert model for chosen direction
    [contrastFromSim,stimulusFromSim] = tfeQCMInvertDirection(params1,theDirection,maxResponse/maxResponseFactor);
    [contrastFromFit,stimulusFromFit] = tfeQCMInvertDirection(paramsQCMFit,theDirection,maxResponse/maxResponseFactor);

    % Plot simulated CRF and inverted points
    figure; hold on
    plot(maxContrast*contrastCoding,directionResponses,'ro','MarkerFaceColor','r','MarkerSize',12);
    plot(contrastFromSim,maxResponse/maxResponseFactor,'bo','MarkerFaceColor','b','MarkerSize',8);
    plot(contrastFromFit,maxResponse/maxResponseFactor,'gx','MarkerSize',14); 
    xlabel('Contrast'); ylabel('Response');
    
    % Plot an isoresponse contour of the simualted and fit model
    nTheta = 100;
    directions = UnitCircleGenerate(nTheta);
    [contrasts1,stimuli1] = tfeQCMInvertDirection(params1,directions,params1.crfAmp/3);
    figure; hold on
    plot(stimuli1(1,:),stimuli1(2,:),'r','LineWidth',3);
    [contrastsFit,stimuliFit] = tfeQCMInvertDirection(paramsQCMFit,directions,params1.crfAmp/3);
    plot(stimuliFit(1,:),stimuliFit(2,:),'b','LineWidth',2);
end

