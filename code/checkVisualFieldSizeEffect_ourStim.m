%% The nominal stimuli from our exteriment

% Load measuremets from session.
dataPath = getpref('LFContrastAnalysis','projectPath');
modulationPath = fullfile(dataPath,'MRContrastResponseFunction','DirectionValidationFiles','KAS25','2018-10-13','postExpValidations.mat');
load(modulationPath)

% Get the SDPs
bkgrdSPD = background.SPDdifferentialDesired(:,1);

% contrast scalars
LminusM_ContrastScalar = .45;
LplusM_ContrastScalar = .491;

% Create the L-M SPD for the postitive arm of the modulation 
directionNum = 1;
LminusM_SPD_pos = bkgrdSPD+(LminusM_ContrastScalar*directedDirection{directionNum}.SPDdifferentialDesired(:,1));

% Create the L+M SPD for the postitive arm of the modulation 
directionNum = 3;
LplusM_SPD_pos = bkgrdSPD+(LplusM_ContrastScalar*directedDirection{directionNum}.SPDdifferentialDesired(:,1));

% Make the Cone fundamentals
S = directedDirection{directionNum}.calibration.describe.S;
wavelengthAxis = SToWls(S);
[~,tempFundamentals] = ComputeCIEConeFundamentals(S,2,30,3);
T_energy = EnergyToQuanta(S,tempFundamentals')';
coneFundamentals = T_energy./max(T_energy')';

% Compute the Contrast
backgroundExcitations = coneFundamentals * bkgrdSPD;
stimExcitations = coneFundamentals * LminusM_SPD_pos;
LminusM_contrast = (stimExcitations-backgroundExcitations)./backgroundExcitations;

% Compute the Contrast
stimExcitations = coneFundamentals * LplusM_SPD_pos;
LplusM_contrast = (stimExcitations-backgroundExcitations)./backgroundExcitations;


% Visual Feild Size Spacing
fieldSizes = 1:20;

[LminusM_QcmStim] = computeContrastStimWithCIE(S,fieldSizes,bkgrdSPD',LminusM_SPD_pos');
[LplusM_QcmStim] = computeContrastStimWithCIE(S,fieldSizes,bkgrdSPD',LplusM_SPD_pos');

% Run the QCM for KAS25
analysisParams = getSubjectParams('KAS25');
fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);

% Use paramters from the QCM fit to KAS25 measurement set 1
params.Qvec        = [0.1661 46.1730];
params.crfAmp      = 2.9998;
params.crfExponent = 1.4826;
params.crfSemi     = 0.2817;
params.expFalloff  = 5.0501; % this does nothing and should be removed
params.crfOffset   = -0.1737;

% Create the Stimulus Struct for L-M
stimulusStruct_LminusM.values = LminusM_QcmStim;
stimulusStruct_LminusM.timebase = 1:length(fieldSizes);

% Create the Stimulus Struct for L-M
stimulusStruct_LplusM.values = LplusM_QcmStim;
stimulusStruct_LplusM.timebase = 1:length(fieldSizes);

% Get the HRF
kernel = [];

% Compute the L-M QCM response
LminusM_Response = fitOBJ.computeResponse(params,stimulusStruct_LminusM,kernel);

% Compute the L-M QCM response
LplusM_Response = fitOBJ.computeResponse(params,stimulusStruct_LplusM,kernel);

% Compute the ratio of the L
LMratio = LminusM_Response.values ./ LplusM_Response.values;

%% Plot some stuff
figure;
subplot(1,3,1)
plot(LminusM_Response.values,'k')
xlabel('Eccentricity')
ylabel('QCM Repsonse')
ylim([0 0.6])
title(sprintf('L-M Response %s%% Contrast',num2str(100*round(norm(LminusM_contrast),2))))
subplot(1,3,2)
plot(LplusM_Response.values,'k')
xlabel('Eccentricity')
ylabel('QCM Repsonse')
ylim([0 0.6])
title(sprintf('L+M Response %s%% Contrast',num2str(100*round(norm(LplusM_contrast),2))))
subplot(1,3,3)
plot(LMratio,'k')
xlabel('Eccentricity')
ylabel('Ratio')
title('L-M/L+M Ratio')

function [QcmStim] = computeContrastStimWithCIE(S,fieldSizes,backgroundSPD,stimulusSPD)

% Loop over visual field sizes
for ii = 1:length(fieldSizes)
    
    % Generate cone fundamentals for different visual field size
    [~,T_quantal] = ComputeCIEConeFundamentals(S,fieldSizes(ii),30,3);
    T_energy = EnergyToQuanta(S,T_quantal')';
    cieConeFund = T_energy./max(T_energy')';
    
    % Compute the background activations
    bkgrdActivations = cieConeFund * backgroundSPD';
    
    % compute the stimulus
    stimActivations = cieConeFund * stimulusSPD';

    % Compute L-M Contrast
    contrast(:,ii) = (stimActivations - bkgrdActivations) ./bkgrdActivations;

end

% Turn the contrasts from above loop into inputs to the QCM
% L-M direction
% Get the angle and contrast
[theta, contrast] = cart2pol(contrast(1,:),contrast(2,:));

% X and y component for a unit vector in the direction of theta
[LconesCont, MConesCont] = pol2cart(theta, 1);

% package for the QCM stim
QcmStim = [LconesCont;MConesCont;contrast];

end