function [cPoints, percentCensored] = findCensoredPoints(motionEstimates,varargin)
% Returns the time points that need to be censored for a given set of motion
% estimates.
%
% Syntax:
%    [cPoints, numCPoints] = censorTimePoints(motionEstimates,varargin)
%
% Description:
%   Calculates the framewise dispalcemnt from the input motion esimates and
%   and applies a threshold to find the time points that need to be
%   censored within a given functional time series.
%
% Inputs:
%   motionEstimates       - Matrix (time x directions) of motion estimates.
%                           direction convention is x y z yaw pitch roll.
%
% Outputs:
%   cPoints               - Vector of frames that need to be censored
%   percentCensored       - The percent of time frames that need to be
%                           censored
%
% Optional key/value pairs:
%   plotMotion            - Flag for plotting the framewise displacement
%                           with markers indicating the timepoints to be
%                           censored. (logical)
%   distMetric            - Distance metric used to compute framewise
%                           displacement. Either 'l1' or 'l2'. (string)
%   sphereSize            - Size (in mm) of the radius used to convert
%                           rotation estimates (in deg or rad) into mm.
%                           Default 50 mm (from Power et al 2012).(numeric)
%   threshold             - Thresholded (in mm) used to determine if a time
%                           point is censored. Default 0.50 mm (Power et
%                           al. 2014) (numeric)
%   rotUnits              - Units of the rotation parameters from the
%                           motion etimates. Either 'deg' or 'rad'.
%                           fMRIprep is in radians and HCP is in degrees.
%                           Default is degrees. (string)
%   addBuffer             - Number of frames to censor before and after
%                           censor point. Input is in the form of [A B]
%                           where A = number of frames prior to censor and
%                           B =  the number of frames after to censor.
%                           Default [0 0]. 
% MAB 08/10/2019 -- wrote it.

%% Input Parser
p = inputParser; p.KeepUnmatched = false;
p.addRequired('motionEstimates', @ismatrix);
p.addParameter('plotMotion',false, @islogical);
p.addParameter('distMetric','l2', @isstr);
p.addParameter('sphereSize',50.0, @isnumeric);
p.addParameter('threshold',0.5, @isnumeric);
p.addParameter('rotUnits','deg', @isstr);
p.addParameter('addBuffer',[], @isvector);
p.parse(motionEstimates, varargin{:})

%% Unpack the parser
radius = p.Results.sphereSize;
threshold = p.Results.threshold;

%% Convert rotation units to mm
rotParams = motionEstimates(:,4:6);
if strcmp(p.Results.rotUnits,'deg')
    rotParamsMM=(2*radius*pi/360)*rotParams;
elseif strcmp(p.Results.rotUnits, 'rad')
    rotParamsMM=radius*temp;
else
    error('Rotation units not recognized');
end
motionEstimates(:,4:6)=rotParamsMM;

%% Go from absolute to relavite motion estimates
relativeMotion=diff(motionEstimates);
relativeMotion=[zeros(1,size(relativeMotion,2)); relativeMotion];

%% Calculate framewise displacement
if strcmp(p.Results.distMetric,'l1')
    fwd = vecnorm(relativeMotion,1,2);
elseif strcmp(p.Results.distMetric,'l2')
    fwd = vecnorm(relativeMotion,2,2);
else
    error('Distance metric not recognized');
end

%% Apply threshold
cPoints    = find(fwd > threshold);
numCPoints = length(cPoints);
percentCensored = numCPoints./length(fwd);

%% Buffer censor points
if ~isempty(p.Results.add


%% Plot
if p.Results.plotMotion
yVals = fwd(cPoints);
figure; hold
fig1 = plot(fwd,'k');
fig1 = plot(cPoints,yVals,'b*');
fig1 = line([0,length(fwd)],[threshold,threshold],'Color','red','LineStyle','--')
legend('Framewise Displacement', 'Censored Points', 'Threshold');

set(fig1, 'LineWidth', 1);

hTitle  = title ('Framewise Displacement w/ Censroed Points');
hXLabel = xlabel('Frames'  );
hYLabel = ylabel('Framewise Displacment (mm)');


set([hTitle, hXLabel, hYLabel],'FontName', 'Helvetica');
set([hXLabel, hYLabel,],'FontSize', 14);
set( hTitle, 'FontSize', 14,'FontWeight' , 'bold');

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:.25:max(fwd)+0.1*max(fwd) , ...
  'LineWidth'   , 2         , ...
  'ActivePositionProperty', 'OuterPosition');

set(gcf, 'Color', 'white' );
end