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

% Plot it
figure; hold on 
plot(wavelengthAxis,backgroundSPD,'k','LineWidth',2)
plot(wavelengthAxis,LminusM_SPDs,'r--','LineWidth',2)
plot(wavelengthAxis,LplusM_SPDs,'b--','LineWidth',2)
legend('background','L-M','L+M')
xlabel('Wavelength (nm)')
ylabel('Power')

