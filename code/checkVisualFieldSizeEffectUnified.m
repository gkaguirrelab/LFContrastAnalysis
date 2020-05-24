% Evaluate how contrasts and predicted MR response vary with eccentricity.
% This version is for our stimuli.

% We'll want to know this at some point
%     theModulation.ConeDirectedDirections{directionNum}.describe.validation
% 
%     ans = 
% 
%     1×10 struct array with fields:
% 
%     background
%     time
%     differentialPrimaryValues
%     measuredPrimaryValues
%     receptors
%     excitationDesired
%     excitationActual
%     contrastDesired
%     contrastActual
%     luminanceDesired
%     luminanceActual
%     label
%     SPDbackground
%     SPDcombined
%     SPDdifferential

% History:
%   05/xx/20  mab  Wrote it.
%   05/18/20  dhb  Adding comments, review

% Clear and close
clear; close all;

% Observer age for calcs
observerAge = 25;
pupilSize = 6;

% Which stimuli
stimulusType = 'monitor';
stimulusType = 'onelight';
switch (stimulusType)
    case 'monitor'
        % Load calibration data and get latest calibration
        load('PTB3TestCal.mat','cals');
        [calStructOBJ, ~] = ObjectToHandleCalOrCalStruct(cals{end});
        
        % Extract the spectral power distributions of the display's RGB primaries
        displaySPDs = (calStructOBJ.get('P_device'));
        
        % Get the cone fundamentals
        S = calStructOBJ.get('S');
        wavelengthAxis = SToWls(S);
        coneFundamentals2deg =  GetHumanPhotoreceptorSS(S, [], 2, observerAge, pupilSize);
        
        % Speficy primary values for background
        bkgrdPrimaries = [0.5 0.5 0.5]';
        
        % Specify desired stimulus contrasts
        LminusM_Modulation = [0.065 -0.065 0]';
        LplusM_Modulation = [0.41 0.41 0]';
        
        % Compute cone excitations for the background
        bkgrdSPD = displaySPDs*bkgrdPrimaries;
        bkgrdExcitations = coneFundamentals2deg * bkgrdSPD;
        
        % Calculate the M matrix
        M = coneFundamentals2deg*displaySPDs;
        
        % Calulate the stimulus excitations
        LminusM_Excitations =  bkgrdExcitations .* (1 + LminusM_Modulation);
        LplusM_Excitations =  bkgrdExcitations .* (1 + LplusM_Modulation);
        
        % Get the stimulus primaries
        LminusM_Primaries = M\LminusM_Excitations;
        LplusM_Primaries = M\LplusM_Excitations;
        
        % Get the stimulus SPDs
        LminusM_SPD_pos = displaySPDs*LminusM_Primaries;
        LplusM_SPD_pos = displaySPDs*LplusM_Primaries;
        
    case 'onelight' 
        ORIGINAL = false;
        if (ORIGINAL)
            % Load the nominal stimuli from our exteriment KAS25 Original session 1
            dataPath = getpref('LFContrastAnalysis','projectPath');
            modulationPath = fullfile(dataPath,'MRContrastResponseFunction','DirectionValidationFiles','KAS25','2018-10-13','postExpValidations.mat');
            theModulation = load(modulationPath);
        else
            % Load test modulation generated more recently
            modulationPath = fullfile(getpref('LFContrastAnalysis','dataDevPath'),'DirectionCorrectedPrimaries','test_mb_good','2020-05-18','correctedPrimaries.mat');
            theModulation = load(modulationPath);
            theModulation.background = theModulation.ConeDirectedBackground;
            theModulation.directedDirection = theModulation.ConeDirectedDirections;
        end
        
        % Get the SDPs of the background
        bkgrdSPD = theModulation.background.SPDdifferentialDesired(:,1);
        bkgrdSPD = bkgrdSPD + theModulation.background.calibration.computed.pr650MeanDark(:,1);
        % bkgrdSPD = theModulation.ConeDirectedBackground.SPDdifferentialDesired(:,1);
        % bkgrdSPD = bkgrdSPD + theModulation.ConeDirectedBackground.calibration.computed.pr650MeanDark(:,1);
        
        % Set the contrast scalars. this is relative to the max contrast of stimuli used in
        % the experiment.  These are set by hand to put responses for two
        % directions in a similar regime.
        % LminusM_ContrastScalar = .45;
        % LplusM_ContrastScalar = .491;
        LminusM_ContrastScalar = 1;
        LplusM_ContrastScalar = 1;
        
        % Get the L-M SPD for the postitive arm of the modulation
        directionNum = 1;
        LminusM_SPD_pos = bkgrdSPD+(LminusM_ContrastScalar*theModulation.directedDirection{directionNum}.SPDdifferentialDesired(:,1));
        %LminusM_SPD_pos = bkgrdSPD+(LminusM_ContrastScalar*theModulation.ConeDirectedDirections{directionNum}.SPDdifferentialDesired(:,1));
        
        % Get the L+M SPD for the postitive arm of the modulation
        directionNum = 3;
        LplusM_SPD_pos = bkgrdSPD+(LplusM_ContrastScalar*theModulation.directedDirection{directionNum}.SPDdifferentialDesired(:,1));
        %LplusM_SPD_pos = bkgrdSPD+(LminusM_ContrastScalar*theModulation.ConeDirectedDirections{directionNum}.SPDdifferentialDesired(:,1));
        
        % Make the cone fundamentals for 2 degrees
        %S = theModulation.directedDirection{directionNum}.calibration.describe.S;
        S = theModulation.directedDirection{directionNum}.calibration.describe.S;

        if (observerAge ~= theModulation.directedDirection{directionNum}.describe.observerAge)
            error('Inconsistent observer age');
        end
        
        % Cone fundamentals
        wavelengthAxis = SToWls(S);
        coneFundamentals2deg =  GetHumanPhotoreceptorSS(S, [], 2, observerAge, pupilSize);
    otherwise
        error('Unknown stimulus type specified');
end

% Compute the L-M Contrast
bkgrdExcitations = coneFundamentals2deg * bkgrdSPD;
stimExcitations = coneFundamentals2deg * LminusM_SPD_pos;
LminusM_contrast = (stimExcitations-bkgrdExcitations)./bkgrdExcitations;

% Compute the L+M Contrast
stimExcitations = coneFundamentals2deg * LplusM_SPD_pos;
LplusM_contrast = (stimExcitations-bkgrdExcitations)./bkgrdExcitations;

% Visual field size spacing
fieldSizes = 1:20;

% Function to loop over the visual field sizes and compute the contrasts and
% and package into QCM friendly intputs
[LminusM_QcmStim,LminusM_fsContrast] = computeContrastStimWithCIE(S,fieldSizes,observerAge,pupilSize,bkgrdSPD',LminusM_SPD_pos');
[LplusM_QcmStim,LplusM_fsContrast] = computeContrastStimWithCIE(S,fieldSizes,observerAge,pupilSize,bkgrdSPD',LplusM_SPD_pos');

% Plot the L and M contrast changes with visual field size
figure; hold on;
s = scatter(LminusM_fsContrast(1,:),LminusM_fsContrast(2,:),90,'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
s.MarkerFaceAlpha = .4;
s = scatter(LplusM_fsContrast(1,:),LplusM_fsContrast(2,:),90,'MarkerEdgeColor',[.5 0 .5],...
    'MarkerFaceColor',[.7 0 .7],'LineWidth',1.5);
s.MarkerFaceAlpha = .4;
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

%% Report 2 and 15 degree contrasts we get here and compare with validation struct
index = find(fieldSizes == 2);
fprintf('L-M: L, M, S contrasts here, 2 deg:                     %0.4f, %0.4f, %0.4f\n',LminusM_fsContrast(1,index),LminusM_fsContrast(2,index),LminusM_fsContrast(3,index));
directionNum = 1;
contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(1:3,1);
fprintf('L-M: L, M, S contrasts nominal from validation, 2 deg:  %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));

index = find(fieldSizes == 15);
fprintf('L-M: L, M, S contrasts here, 15 deg:                    %0.4f, %0.4f, %0.4f\n',LminusM_fsContrast(1,index),LminusM_fsContrast(2,index),LminusM_fsContrast(3,index));
directionNum = 1;
contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(4:6,1);
fprintf('L-M: L, M, S contrasts nominal from validation, 15 deg: %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));
fprintf('\n')

index = find(fieldSizes == 2);
fprintf('L+M: L, M, S contrasts here, 2 deg:                     %0.4f, %0.4f, %0.4f\n',LplusM_fsContrast(1,index),LplusM_fsContrast(2,index),LplusM_fsContrast(3,index));
directionNum = 3;
contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(1:3,1);
fprintf('L+M: L, M, S contrasts nominal from validation, 15 deg: %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));

index = find(fieldSizes == 15);
fprintf('L+M: L, M, S contrasts here, 15 deg:                    %0.4f, %0.4f, %0.4f\n',LplusM_fsContrast(1,index),LplusM_fsContrast(2,index),LplusM_fsContrast(3,index));
directionNum = 3;
contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(4:6,1);
fprintf('L+M: L, M, S contrasts nominal from validation, 15 deg: %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));

% Sub-function
function [QcmStim,contrasts] = computeContrastStimWithCIE(S,fieldSizes,age,pupilSize,backgroundSPD,stimulusSPD)

% Loop over visual field sizes
for ii = 1:length(fieldSizes)
    
    % Generate cone fundamentals for different visual field size
    cieConeFund =  GetHumanPhotoreceptorSS(S, [], fieldSizes(ii), age,pupilSize);
    
    % Compute the background activations
    bkgrdActivations = cieConeFund * backgroundSPD';
    
    % Compute the stimulus
    stimActivations = cieConeFund * stimulusSPD';
    
    % Compute L-M Contrast
    contrasts(:,ii) = (stimActivations - bkgrdActivations) ./ bkgrdActivations;
    
end

% Turn the contrasts from above loop into inputs to the QCM
% L-M direction
[theta, contrastVals] = cart2pol(contrasts(1,:),contrasts(2,:));
[LconesCont, MConesCont] = pol2cart(theta, 1);

% package for the QCM stim
QcmStim = [LconesCont;MConesCont;contrastVals];

end