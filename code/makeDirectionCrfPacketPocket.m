function directionCrfMeanPacketPocket = makeDirectionCrfPacketPocket(analysisParams,iampParams, varargin);
% Takes packets procuded by fit_IAMP and replaces the stimulus with 
% direction and contrast for contrast response function fits.
%
% Syntax:
%    directionTimeCoursePacketPocket = makeDirectionTimeCoursePacketPocket(iampTimeCoursePacketPocket)
%
% Description:
%    Takes packets procuded by fit_IAMP and replaces the existing stimulus
%    that is coded for for the IAMP model with the direction and contrast
%    stimulus need for both the Naka-Rushton and QCM models. The contrast
%    and directions information lives in the metaData subfield of the IAMP
%    packets
%
% Inputs:
%    analysisParams              - Struct of important information for the
%                                  analysis 
%
%    avgIampParams                - Cell array of IAMP packets from fit_IAMP 
%
% Outputs:
%    directionCrfMeanPacketPocket - Cell array of packets for directions
%                                   based models such as QCM and Naka-Ruston 
%
% Optional key/value pairs:
%    none

% MAB 01/06/19

p = inputParser;
p.addRequired('analysisParams',@isstruct);
p.addRequired('iampParams',@isstruct);
p.parse(analysisParams, iampParams,varargin{:});

% Construct the stimulus 
stim = kron(analysisParams.directionCoding(1:analysisParams.theDimension,:).*analysisParams.maxContrastPerDir,analysisParams.contrastCoding);
stim = [stim, [0;0]];
[qcmDirStimDirections,qcmDirStimContrasts] = tfeQCMStimuliToDirectionsContrasts(stim,'precision',4);

% Make the packet - Stimulus
directionCrfMeanPacketPocket.stimulus.values   = [qcmDirStimDirections; qcmDirStimContrasts];
directionCrfMeanPacketPocket.stimulus.timebase = 1:size(directionCrfMeanPacketPocket.stimulus.values,2);

% Make the packet - Response 
directionCrfMeanPacketPocket.response.values = iampParams.paramMainMatrix;
directionCrfMeanPacketPocket.response.timbase = directionCrfMeanPacketPocket.stimulus.timebase;

% Add an empty kernel and metadata 
directionCrfMeanPacketPocket.kernel            = [];
directionCrfMeanPacketPocket.metaData          = [];

end


