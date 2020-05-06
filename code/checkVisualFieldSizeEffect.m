% Set paths to the cal data and con fundamentals
[~, userID] = system('whoami');
userID = strtrim(userID);
pathToPTB = fullfile('/Users',userID,'Documents','MATLAB','toolboxes','Psychtoolbox-3','Psychtoolbox');
pathToCalData = fullfile(pathToPTB ,'PsychCalDemoData','PTB3TestCal.mat');
pathToConeFund = fullfile(pathToPTB,'PsychColorimetricData','PsychColorimetricMatFiles','T_cones_ss2.mat');

% Load the Smith-Pokorny 2 deg cone fundamentals
load(pathToConeFund,'S_cones_ss2','T_cones_ss2');

% Load Cal Data
load(pathToCalData,'cals');
% Get the latest calibration
[calStructOBJ, ~] = ObjectToHandleCalOrCalStruct(cals{end});

% Extract the spectral power distributions of the display's RGB primaries
displaySPDs = (calStructOBJ.get('P_device'))';

% Spline the Smith-Pokorny 2 deg cone fundamentals to match the wavelengthAxis
S = calStructOBJ.get('S');
wavelengthAxis = SToWls(S);
coneFundamentals = SplineCmf(S_cones_ss2, T_cones_ss2, WlsToS(wavelengthAxis));

% Speficy primary values for background
backgroundPrimaries = [0.5 0.5 0.5]';

% Speficy stimulus contrast ~10%
LminusM_Modulation = [0.07 -0.07 0]';
LplusM_Modulation = [0.07 0.07 0]';

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

% % Plot it
% figure; hold on
% plot(wavelengthAxis,backgroundSPD,'k','LineWidth',2)
% plot(wavelengthAxis,LminusM_SPDs,'r--','LineWidth',2)
% plot(wavelengthAxis,LplusM_SPDs,'b--','LineWidth',2)
% legend('background','L-M','L+M')
% xlabel('Wavelength (nm)')
% ylabel('Power')

% Visual Feild Size Spacing
fieldSizes = 1:20;

% Loop over visual field sizes
for ii = 1:length(fieldSizes)
    % generate cone fundamentals for different visual field size
    cieConeFund = ComputeCIEConeFundamentals(S,fieldSizes(ii),30,3);
    
    bkgrd = cieConeFund * backgroundSPD';
    A = cieConeFund * LminusM_SPDs';
    B = cieConeFund * LplusM_SPDs';
    
    LminusM_Contrast(:,ii) = (A - bkgrd) ./bkgrd;
    LplusM_Contrast(:,ii)  = (B - bkgrd) ./bkgrd;
    
end

% Turn the Contrast from above loop into inputs to the QCM
% L-M direction
% Get the angle and contrast
[theta, contrast] = cart2pol(LminusM_Contrast(1,:),LminusM_Contrast(2,:));
% X and y component for a unit vector in the direction of theta
[LconesCont, MConesCont] = pol2cart(theta, 1);
% package for the QCM stim
LminusM_QcmStim = [LconesCont;MConesCont;contrast];

% L+M direction
% Get the angle and contrast
[theta, contrast] = cart2pol(LplusM_Contrast(1,:),LplusM_Contrast(2,:));
% X and y component for a unit vector in the direction of theta
[LconesCont, MConesCont] = pol2cart(theta, 1);
% package for the QCM stim
LplusM_QcmStim = [LconesCont;MConesCont;contrast];

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
[analysisParams] = loadHRF(analysisParams);
kernel = analysisParams.HRF;

% Compute the L-M QCM response
LminusM_Response = fitOBJ.computeResponse(params,stimulusStruct_LminusM,kernel);

% Compute the L-M QCM response
LplusM_Response = fitOBJ.computeResponse(params,stimulusStruct_LplusM,kernel);

% Compute the ratio of the L
LMratio = LminusM_Response.values ./ LplusM_Response.values;
