% Load the nominal stimuli from our exteriment KAS25 Original session 1
% dataPath = getpref('LFContrastAnalysis','projectPath');
% modulationPath = fullfile(dataPath,'MRContrastResponseFunction','DirectionValidationFiles','KAS25','2018-10-13','postExpValidations.mat');
modulationPath = fullfile('/Users','michael','labDropbox','MELA_datadev','Experiments','OLApproach_TrialSequenceMR',...
         'MRCRF','DirectionCorrectedPrimaries','test_mb_good','2020-05-18','correctedPrimaries.mat');
load(modulationPath)

% Get the SDPs of the background
bkgrdSPD = ConeDirectedBackground.SPDdifferentialDesired(:,1);
bkgrdSPD = bkgrdSPD + ConeDirectedBackground.calibration.computed.pr650MeanDark(:,1);
% Set the contrast scalars. this is relative to the max contrast of stimuli used in
% the experiment 
LminusM_ContrastScalar = 1;
LplusM_ContrastScalar = 1;

% Create the L-M SPD for the postitive arm of the modulation
directionNum = 1;
LminusM_SPD_pos = bkgrdSPD+(LminusM_ContrastScalar*ConeDirectedDirections{directionNum}.SPDdifferentialDesired(:,1));

% Create the L+M SPD for the postitive arm of the modulation
directionNum = 3;
LplusM_SPD_pos = bkgrdSPD+(LplusM_ContrastScalar*ConeDirectedDirections{directionNum}.SPDdifferentialDesired(:,1));

% Make the Cone fundamentals for 2 degrees
observerAgeInYears = 25;
S = ConeDirectedDirections{directionNum}.calibration.describe.S;
wavelengthAxis = SToWls(S);
lambdaMaxShift = [];
pupilDiameterMm = 8;
fieldSizeDegrees = [2 2 2];
photoreceptorClasses = {'LConeTabulatedAbsorbance','MConeTabulatedAbsorbance','SConeTabulatedAbsorbance'};
fractionBleached = OLEstimateConePhotopigmentFractionBleached(S,bkgrdSPD,pupilDiameterMm,fieldSizeDegrees,observerAgeInYears,photoreceptorClasses);
coneFundamentals = GetHumanPhotoreceptorSS(S,photoreceptorClasses,fieldSizeDegrees,observerAgeInYears,pupilDiameterMm,lambdaMaxShift,fractionBleached);
%coneFundamentals =  GetHumanPhotoreceptorSS(S, [], 2, 25);

% Compute the L-M Contrast
backgroundExcitations = coneFundamentals * bkgrdSPD;
stimExcitations = coneFundamentals * LminusM_SPD_pos;
LminusM_contrast = (stimExcitations-backgroundExcitations)./backgroundExcitations;

% Compute the L+M Contrast
stimExcitations = coneFundamentals * LplusM_SPD_pos;
LplusM_contrast = (stimExcitations-backgroundExcitations)./backgroundExcitations;


% Visual Feild Size Spacing
fieldSizes = 1:20;

% function to loop over the visual field sizes and compute the contrasts and
% and package into QCM friendly intputs 
[LminusM_QcmStim,LminusM_fsContrast] = computeContrastStimWithCIE(S,fieldSizes,bkgrdSPD',LminusM_SPD_pos');
[LplusM_QcmStim,LplusM_fsContrast] = computeContrastStimWithCIE(S,fieldSizes,bkgrdSPD',LplusM_SPD_pos');

% plot the L and M contrast changes with visual field size
figure; hold on;
s = scatter(LminusM_fsContrast(1,:),LminusM_fsContrast(2,:),90,'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
s.MarkerFaceAlpha = .4;
s = scatter(LplusM_fsContrast(1,:),LplusM_fsContrast(2,:),90,'MarkerEdgeColor',[.5 0 .5],...
    'MarkerFaceColor',[.7 0 .7],'LineWidth',1.5);
s.MarkerFaceAlpha = .4
line([-.1 .3], [0 0], 'Color', [.3 .3 .3], 'LineStyle', ':','LineWidth', 2);
line([0 0], [-.1 .3], 'Color', [.3 .3 .3], 'LineStyle', ':','LineWidth', 2);
axis square
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
legend('L minus M','L plus M')

%% Run the QCM for KAS25
% Load exp params
analysisParams = getSubjectParams('KAS25');

% Create QCM object
fitOBJ = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);

% Use parameters from the QCM fit to KAS25 measurement set 1
params.Qvec        = [0.1661 46.1730];
params.crfAmp      = 2.9998;
params.crfExponent = 1.4826;
params.crfSemi     = 0.2817;
params.expFalloff  = 5.0501; % this does nothing and should be removed for fit
params.crfOffset   = -0.1737;

% Create the Packet for L-M
stimulusStruct_LminusM.values = LminusM_QcmStim;
stimulusStruct_LminusM.timebase = 1:length(fieldSizes);

% Create the Packet for L+M
stimulusStruct_LplusM.values = LplusM_QcmStim;
stimulusStruct_LplusM.timebase = 1:length(fieldSizes);

% Get the HRF
kernel = [];

% Compute the L-M QCM response
LminusM_Response = fitOBJ.computeResponse(params,stimulusStruct_LminusM,kernel);

% Compute the L+M QCM response
LplusM_Response = fitOBJ.computeResponse(params,stimulusStruct_LplusM,kernel);

% Compute the ratio of the L-M/L+M
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

% Sub-function 
function [QcmStim,contrasts] = computeContrastStimWithCIE(S,fieldSizes,backgroundSPD,stimulusSPD)

% Loop over visual field sizes
for ii = 1:length(fieldSizes)
    
    % Generate cone fundamentals for different visual field size
    cieConeFund =  GetHumanPhotoreceptorSS(S, [], fieldSizes(ii), 25);
    
    
    % Compute the background activations
    bkgrdActivations = cieConeFund * backgroundSPD';
    
    % compute the stimulus
    stimActivations = cieConeFund * stimulusSPD';
    
    % Compute L-M Contrast
    contrasts(:,ii) = (stimActivations - bkgrdActivations) ./bkgrdActivations;
    
end

% Turn the contrasts from above loop into inputs to the QCM
% L-M direction
% Get the angle and contrast
[theta, contrastVals] = cart2pol(contrasts(1,:),contrasts(2,:));

% X and y component for a unit vector in the direction of theta
[LconesCont, MConesCont] = pol2cart(theta, 1);

% package for the QCM stim
QcmStim = [LconesCont;MConesCont;contrastVals];

end