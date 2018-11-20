%--------------------------------------------------------------------------
%% Set up params:
%--------------------------------------------------------------------------
% Naka Ruston Params
Rmax   = .9;
sigma  = .7;
n      = 2.1;
offset = 0;

% Ellipse Params
minorAxis = .3;
rotdeg    = 45;

% Param Vec
ellParams = [1, minorAxis, rotdeg];


%--------------------------------------------------------------------------
%% Generate Stimulus:
%--------------------------------------------------------------------------
% Random stimuli bounded between 1 and -1
numStim = 300;
stim = (2*rand(2,numStim) -1);

%--------------------------------------------------------------------------
%% Generate Response:
%--------------------------------------------------------------------------
% Set up scale matrix
S = diag(ellParams(1:2));
% Set up rotation matrix
V = deg2rotm(ellParams(3))';  

% Get the Q matrix that takes stim and get a radius. 
A = S*V';
Q = A'*A;

% Get the raduis
radius =  diag(sqrt(stim'*Q*stim));

% Get the neural response values
R  = nakaRushton(radius,sigma,n,Rmax, offset);

%--------------------------------------------------------------------------
%%  Use the tfeQCM to fit the stim/resp:
%--------------------------------------------------------------------------
% Get the tfeQCM object
temporalFitQCM = tfeQCM('verbosity','none','dimension',2);

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
[paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket);
fprintf('Model parameter from fits:\n');
temporalFitQCM.paramPrint(paramsQCMFit)

%--------------------------------------------------------------------------
%%  Visualize the effects of the fit
%--------------------------------------------------------------------------

qcmParams =[1 paramsQCMFit.Qvec];
S_qcm = diag(qcmParams(1:2));
% Set up rotation matrix
V_qcm = deg2rotm(qcmParams(3))';  

% Get the Q matrix that takes stim and get a radius. 
A_qcm = S_qcm*V_qcm';
Q_qcm = A_qcm'*A_qcm;

% Get the raduis
radiusQcm =  diag(sqrt(stim'*Q_qcm*stim));

% Get the neural response values
Rqcm  = nakaRushton(radiusQcm,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent,paramsQCMFit.crfAmp, paramsQCMFit.offset);

%--------------------------------------------------------------------------
%%  Invert the model
%--------------------------------------------------------------------------
thresh = 0.3;
eqContrast = InvertNakaRushton([paramsQCMFit.crfAmp,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent],thresh);
circlePoints = eqContrast*UnitCircleGenerate(numStim);
[~,Ainv,Q] = EllipsoidMatricesGenerate([1 paramsQCMFit.Qvec],'dimension',2);
ellipsePoints = Ainv*circlePoints;
checkThresh = ComputeNakaRushton([paramsQCMFit.crfAmp,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent],diag(sqrt(ellipsePoints'*Q*ellipsePoints)));
if (any(abs(checkThresh-thresh) > 1e-10))
    error('Did not invert QCM model correctly');
end

%--------------------------------------------------------------------------
%%  Plot
%--------------------------------------------------------------------------
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





