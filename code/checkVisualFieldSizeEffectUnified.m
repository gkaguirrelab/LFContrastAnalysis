% Evaluate how contrasts and predicted MR response vary with eccentricity.
% This version is for our stimuli.


% This info might be useful, from saveParameters function in experimental
% code.
%
% %set the contrast vector angle in the LM plane.
% %  0 = L, 90 = M, 45 = L+M, -45 = L-M;   [-22.5, 22.5, 67.5, 112.5];
% protocolParams.LMVectorAngles =[-45, 0, 45, 90];
% protocolParams.contrastLevels =  trialTypeParams.contrastLevels;
% % order below matches the contrast angle vector above.
% % Max values for Box D in 10/18 are [-22.5,..,112.5] => [0.085,0.20,0.40,0.13] and [-45,...,90] => [.12, .14, .60 , .22];

% History:
%   05/xx/20  mab  Wrote it.
%   05/18/20  dhb  Adding comments, review

% Clear and close
clear; close all;

% Observer age for calcs
observerAge = 25;
pupilSize = 6;

% Set the contrast scalars. this is relative to the max contrast of stimuli used in
% the experiment.
LminusM_ContrastScalar = 1;
LplusM_ContrastScalar = 1;

% Which stimulus type?
% This can be
%    stimulusType = 'monitor';
%    stimulusType = 'onelight';
stimulusType = 'onelight';
switch (stimulusType)
    case 'monitor'
        % Load some monitor calibration data and get latest calibration
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
        LminusM_Modulation = [0.085 -0.085 0]';
        LplusM_Modulation = [0.40 0.40 0]';
        
        % Compute cone excitations for the background
        bkgrdSPD = displaySPDs*bkgrdPrimaries;
        bkgrdExcitations = coneFundamentals2deg * bkgrdSPD;
        
        % Calculate the M matrix
        M = coneFundamentals2deg*displaySPDs;
        
        % Calulate the stimulus excitations
        LminusM_Excitations =  bkgrdExcitations .* (1 + LminusM_ContrastScalar*LminusM_Modulation);
        LplusM_Excitations =  bkgrdExcitations .* (1 + LplusM_ContrastScalar*LplusM_Modulation);
        
        % Get the stimulus primaries
        LminusM_Primaries = M\LminusM_Excitations;
        LplusM_Primaries = M\LplusM_Excitations;
        
        % Get the stimulus SPDs
        LminusM_SPD_pos = displaySPDs*LminusM_Primaries;
        LplusM_SPD_pos = displaySPDs*LplusM_Primaries;
        
    case 'onelight'
        % ORIGINAL can be false or true. If false, read a file as in the
        % experiment.  If true, read a file that we created more recently
        % as a test.  The test does not have actual validation data.
        ORIGINAL = true;
        if (ORIGINAL)
            % Load the nominal stimuli from our exteriment KAS25,
            % measurement set 1. Session 1, 2018-10-13, Session 2, 2018-10-20
            dataPath = getpref('LFContrastAnalysis','projectPath');
            modulationPath = fullfile(dataPath,'MRContrastResponseFunction','DirectionValidationFiles','KAS25','2018-13-20','postExpValidations.mat');
            theModulation = load(modulationPath);
        else
            % Load test modulation generated more recently.  Rename some
            % variables to match what happens when we load an original
            % file.
            modulationPath = fullfile(getpref('LFContrastAnalysis','dataDevPath'),'DirectionCorrectedPrimaries','test_mb_good','2020-05-18','correctedPrimaries.mat');
            theModulation = load(modulationPath);
            theModulation.background = theModulation.ConeDirectedBackground;
            theModulation.directedDirection = theModulation.ConeDirectedDirections;
        end
        
        % Get the SDPs of the background
        bkgrdSPD = theModulation.background.SPDdifferentialDesired(:,1);
        bkgrdSPD = bkgrdSPD + theModulation.background.calibration.computed.pr650MeanDark(:,1);
        
        % Get the L-M SPD for the postitive arm of the modulation
        directionNum = 1;
        LminusM_SPD_pos = bkgrdSPD+(LminusM_ContrastScalar*theModulation.directedDirection{directionNum}.SPDdifferentialDesired(:,1));
        
        % Get the L+M SPD for the postitive arm of the modulation
        directionNum = 3;
        LplusM_SPD_pos = bkgrdSPD+(LplusM_ContrastScalar*theModulation.directedDirection{directionNum}.SPDdifferentialDesired(:,1));
        
        % Make the cone fundamentals for 2 degrees
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
xlim([-0.2 0.6]); ylim([-0.2 0.6]);
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

% Create the packet for L-M
stimulusStruct_LminusM.values = LminusM_QcmStim;
stimulusStruct_LminusM.timebase = 1:length(fieldSizes);

% Create the packet for L+M
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

%% Plot predicted response and response ratio
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
ylim([0 15]);
title('L-M/L+M Ratio')

%% Report 2 and 15 degree contrasts we get here and compare with validation struct
index = find(fieldSizes == 2);
fprintf('L-M: L, M, S contrasts here, 2 deg:                     %0.4f, %0.4f, %0.4f\n',LminusM_fsContrast(1,index),LminusM_fsContrast(2,index),LminusM_fsContrast(3,index));
if (strcmp(stimulusType,'onelight'))
    directionNum = 1;
    contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(1:3,1);
    fprintf('L-M: L, M, S contrasts nominal from validation, 2 deg:  %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));
end

index = find(fieldSizes == 15);
fprintf('L-M: L, M, S contrasts here, 15 deg:                    %0.4f, %0.4f, %0.4f\n',LminusM_fsContrast(1,index),LminusM_fsContrast(2,index),LminusM_fsContrast(3,index));
if (strcmp(stimulusType,'onelight'))
    directionNum = 1;
    contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(4:6,1);
    fprintf('L-M: L, M, S contrasts nominal from validation, 15 deg: %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));
end
fprintf('\n')

index = find(fieldSizes == 2);
fprintf('L+M: L, M, S contrasts here, 2 deg:                     %0.4f, %0.4f, %0.4f\n',LplusM_fsContrast(1,index),LplusM_fsContrast(2,index),LplusM_fsContrast(3,index));
if (strcmp(stimulusType,'onelight'))
    directionNum = 3;
    contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(1:3,1);
    fprintf('L+M: L, M, S contrasts nominal from validation, 15 deg: %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));
end

index = find(fieldSizes == 15);
fprintf('L+M: L, M, S contrasts here, 15 deg:                    %0.4f, %0.4f, %0.4f\n',LplusM_fsContrast(1,index),LplusM_fsContrast(2,index),LplusM_fsContrast(3,index));
if (strcmp(stimulusType,'onelight'))
    directionNum = 3;
    contrastValidation = theModulation.directedDirection{directionNum}.describe.validation(1).contrastDesired(4:6,1);
    fprintf('L+M: L, M, S contrasts nominal from validation, 15 deg: %0.4f, %0.4f, %0.4f\n',contrastValidation(1),contrastValidation(2),contrastValidation(3));
end

%% Compute contrasts from validation measurements
for fieldSize = [2,15]
    switch (fieldSize)
        case 2
            fieldSizeRowIndices = 1:3;
        case 15
            fieldSizeRowIndices = 4:6;
        otherwise
            error('Comparison with original requires field size of 2 or 15');
    end
    
    % Go through each validation measurement and get positive and negative
    % contrast.  Does some sanity checking that our calculations now match
    % those done at validation time.
    Tcones = GetHumanPhotoreceptorSS(S, [], fieldSize, observerAge, pupilSize);
    for posOrNeg = 1:2
        for dd = 1:length(theModulation.directedDirection)
            for ii = 1:length(theModulation.directedDirection{1}.describe.validation)
                measuredBackgroundSpd(:,ii,dd) = theModulation.directedDirection{dd}.describe.validation(ii).SPDbackground.measuredSPD;
                measuredModulationSpd(:,ii,dd) = theModulation.directedDirection{dd}.describe.validation(ii).SPDcombined(:,posOrNeg).measuredSPD;
                lmsBackground(:,ii,dd) = Tcones*measuredBackgroundSpd(:,ii,dd);
                lmsModulation(:,ii,dd) = Tcones*measuredModulationSpd(:,ii,dd);
                lmsContrast(:,ii,dd) = ExcitationsToContrast(lmsModulation(:,ii,dd),lmsBackground(:,ii,dd));
                lmsContrastFromBefore(:,ii,dd) = theModulation.directedDirection{dd}.describe.validation(ii).contrastActual(fieldSizeRowIndices,posOrNeg);
                if (posOrNeg == 1)
                    lmsContrastPos(:,ii,dd) = lmsContrast(:,ii,dd);
                else
                    lmsContrastNeg(:,ii,dd) = lmsContrast(:,ii,dd);
                end
            end
        end
        if (max(abs(lmsContrast(:) - lmsContrastFromBefore(:))) > 1e-10)
            error('Not reproducing contrast calculation');
        end
    end
    
    % Get median positive and negative contrast for each direction,
    % and mean of the two.  Negative contrasts have sign flipped so
    % all are expressed using the postive convention.
    for dd = 1:length(theModulation.directedDirection)
        posContrasts = squeeze(lmsContrastPos(:,:,dd));
        medianPosContrasts(:,dd) = median(posContrasts,2);
        negContrasts = squeeze(lmsContrastNeg(:,:,dd));
        medianNegContrasts(:,dd) = -median(negContrasts,2);
        meanContrasts(:,dd) = mean([medianPosContrasts(:,dd), medianNegContrasts(:,dd)],2);
    end
    
    fprintf('\n****\nField size %d\n****\n',fieldSize)
    fprintf('Median validation positive contrasts\n');
    medianPosContrasts
    fprintf('Median validation negative contrasts (sign flipped)\n');
    medianNegContrasts
    fprintf('Mean of validation pos and neg contrasts\n');
    meanContrasts
    fprintf('\n\n')
end


%% Sub-function
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