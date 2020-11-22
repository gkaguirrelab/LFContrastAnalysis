function [equivContrasts] = computeEquivContrast(stimPacket,qcmParams,varargin)
%  Use the paramters of the QCM to convert the stimulus into equivalent
%  contrast.
%
% Syntax:
%   [equivContrasts] = computeEquivContrast(stimPacket,qcmParams,varargin)
%
% Description:
%
%
% Inputs:
%    stimPacket        - A 3xN matrix specifying the stimuli:
%                        Rows 1&2 = the x and y vector components of a
%                        unit vector in the stimulis direction.
%                        Row 3    = the corresponding vector magnitude.
%    qcmParams         - Parameters from the QCM model fit.
%                           Angle, Minor Axis Ratio, NR Exp, NR AMP, NR
%                           Semi, NR Offset
%
% Outputs:
%    equivContrasts    - The equivalent contrasts. Row 1 of output.
%
% Optional key/value pairs:
%    addIampResp       - pair the equivalent contrast with its
%                        corresponding IAMP beta weight. Row 2 of output.
%    addPlotColor      - This appends a custon color scheme for the
%                        LFContrast data set. Rows 3-6;

% MAB 03/24/20

%% Input Parser
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('stimPacket',@isstruct);
p.addRequired('qcmParams',@isstruct);
p.addParameter('addIampResp',true,@islogical);
p.addParameter('addPlotColor',true,@islogical);
p.parse(stimPacket,qcmParams,varargin{:});

% Get the stimuli from the packet
stimuli = stimPacket.stimulus.values;

% Get the stumilus directions
direction = stimuli(1:2,:);

% Get the stumilus magnitude
magitude  = stimuli(3,:);

% Get the stimuls vectors in the LM plane
directionVecs = direction .* magitude;

% Get the dimensions of the stimulus
dimension = size(direction,1);

% Generate the Q matrix to transform the stimuli to eqiv. contrast
[~,~,Q] = EllipsoidMatricesGenerate([1 1./qcmParams.Qvec(1) qcmParams.Qvec(2)]','dimension',dimension);

% Transform the stimuli to eqiv. contrast

if p.Results.addIampResp
    equivContrasts.values =[sqrt(diag(directionVecs'*Q*directionVecs))';...
                            stimPacket.response.values];
else
    equivContrasts.values = sqrt(diag(directionVecs'*Q*directionVecs))';

end

%% THIS SECTION IS ONLY RELEVENT TO PLOTTING DATA FOR LFCONTRAST PAPER
if p.Results.addPlotColor
    equivContrasts.colorVals = [repmat([50,136,189],[5,1]);...
                 repmat([254,224,139],[5,1]);...
                 repmat([171,221,164],[5,1]);...
                 repmat([244,109,67],[5,1]);...
                 repmat([102,194,165],[5,1]);...
                 repmat([230,245,152],[5,1]);...
                 repmat([253,174,97],[5,1]);...
                 repmat([213,62,79],[5,1]);...
                 [0,0,0]]./255;
    equivContrasts.colorMap = [repmat([50,136,189],[32,1]);...
                 repmat([102,194,165],[32,1]);...
                 repmat([171,221,164],[32,1]);...
                 repmat([230,245,152],[32,1]);...
                 repmat([254,224,139],[32,1]);...
                 repmat([253,174,97],[32,1]);...
                 repmat([244,109,67],[32,1]);...
                 repmat([213,62,79],[32,1])]./255;
end

end