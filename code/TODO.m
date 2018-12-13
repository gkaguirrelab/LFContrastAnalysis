%% GOALS
%
% 1) Fit Naka-Rushton function to the MR timecourse data, with NR's
% expressed separately for each color direction in the stimulus.
%
% 1a) Explore various constraints (common amplitude, etc.) of the NR
% fits.
%
% 2) Use the fits to be clever about how we choose starting parameters for
% the QCM, and then see if we can fit that directly to the timecourse on
% each run.

%% Now we are going to fit the timecourse with Naka-Rushton crfs
% 
% Use tfeNakaRushtonDirection to fit Naka-Rushton functions to all
% of the directions.

% Make a packet that represents the color direciton and contrast at each
% time point.
LMSContrastMat = LMSContrastValuesFromParams(expParams,analysisParams.contrastCoding,analysisParams.directionCoding,analysisParams.maxContrastPerDir,totalTime,deltaT);