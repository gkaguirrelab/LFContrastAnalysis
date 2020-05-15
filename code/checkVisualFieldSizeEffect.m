% Load Cal Data
load('PTB3TestCal.mat','cals');
% Get the latest calibration
[calStructOBJ, ~] = ObjectToHandleCalOrCalStruct(cals{end});

% Extract the spectral power distributions of the display's RGB primaries
displaySPDs = (calStructOBJ.get('P_device'))';

% Get the cone fundamentals
S = calStructOBJ.get('S');
wavelengthAxis = SToWls(S);
[~,tempFundamentals] = ComputeCIEConeFundamentals(S,2,30,3);
T_energy = EnergyToQuanta(S,tempFundamentals')';
coneFundamentals = T_energy./max(T_energy')';

% Speficy primary values for background
backgroundPrimaries = [0.5 0.5 0.5]';

% Speficy stimulus contrast
LminusM_Modulation = [0.065 -0.065 0]';
LplusM_Modulation = [0.41 0.41 0]';

% Compute cone excitations for the background
backgroundSPD = backgroundPrimaries(1) * displaySPDs(1,:) + ...
    backgroundPrimaries(2) * displaySPDs(2,:) + ...
    backgroundPrimaries(3) * displaySPDs(3,:);

backgroundExcitations = coneFundamentals * backgroundSPD';

% Calculate the M matrix
M = coneFundamentals *displaySPDs';

% Calulate the stimulus excitations
LminusM_Excitations =  backgroundExcitations .* (1 + LminusM_Modulation);
LplusM_Excitations =  backgroundExcitations .* (1 + LplusM_Modulation);

% Get the stimulus primaries
LminusM_Primaries = M\LminusM_Excitations;
LplusM_Primaries = M\LplusM_Excitations;

% Get the stimulus SPDs
LminusM_SPDs = sum(displaySPDs.*LminusM_Primaries,1);
LplusM_SPDs = sum(displaySPDs.*LplusM_Primaries,1);

% Visual Feild Size Spacing
fieldSizes = 1:20;

[LminusM_QcmStim] = computeContrastStimWithCIE(S,fieldSizes,backgroundSPD,LminusM_SPDs);
[LplusM_QcmStim] = computeContrastStimWithCIE(S,fieldSizes,backgroundSPD,LplusM_SPDs);



% Make the stim
[LminusM_QcmStim] = computeContrastStimWithCIE(ConeDirectedDirections{directionNum}.calibration.describe.S,fieldSizes,bkgrdSPD',LminusM_SPD_pos');
[LplusM_QcmStim] = computeContrastStimWithCIE(ConeDirectedDirections{directionNum}.calibration.describe.S,fieldSizes,bkgrdSPD',LplusM_SPD_pos');

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
title(sprintf('L-M Response %s%% Contrast',num2str(100*round(norm(LminusM_Modulation),2))))
subplot(1,3,2)
plot(LplusM_Response.values,'k')
xlabel('Eccentricity')
ylabel('QCM Repsonse')
ylim([0 0.6])
title(sprintf('L+M Response %s%% Contrast',num2str(100*round(norm(LplusM_Modulation),2))))
subplot(1,3,3)
plot(LMratio,'k')
xlabel('Eccentricity')
ylabel('Ratio')
title('L-M/L+M Ratio')

% Visual Feild Size Spacing
fieldSizes = 1:20;

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