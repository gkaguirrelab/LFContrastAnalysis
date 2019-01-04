function [] = plotQCMtimecourse(analysisParams, packets, varargin)
% Takes in a text file name and retuns a cell of the lines of the text file
%
% Syntax:
%   filesCell = plotQCMtimecourse(inFile)
%
% Description:
%    This function takes in a file name for a trext file and returns a cell
%    that is composed of the lines of the text file. Example of this would
%    be a text file of file names the output is a cell of files names.
%
% Inputs:
%    inFile            - File name of a text file. (string)
%
% Outputs:
%    fileCell          - A cell of the lines of the input text file. (cell)
%
% Optional key/value pairs:
%    none

% MAB 09/09/18
% MAB 01/03/19 -- made more general

% paramsFitIAMP,meanIAMPBetas,fitResponseStructQCM,baseline
% variable i need to add as key vaulue pairs plus fitParams

% NOTE: MB: change so the the input for the mean beta weights is the iamp
% fit with the paramsMainMatrix altered to contian the mean betas (as is
% done below. this allows for one less key/value pair.
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('analysisParams',@isstruct);
p.addRequired('packets',@iscell);
p.addParameter('paramsIAMP',[],iscell);
p.addParameter('paramsMeanBetas',[],@(x)(isvector(x) | iscell(x)));
p.addParameter('paramsQCMDir',[],@isstruct);
p.addParamter('paramsNakaRushton',[],@isstruct);
p.addParamter('nakaRushtonDirections',[],@ismatrix)
p.addParamter('baseline',0,@isnumeric);
p.parse(thePacket,varargin{:});

% Sort varargin
paramsIAMP            = p.Results.paramsIAMP;
paramsMeanBetas       = p.Results.paramsMeanBetas;
paramsQCMDir          = p.Results.paramsQCMDir;
paramsNakaRushton     = p.Results.paramsNakaRushton;
nakaRushtonDirections = p.Resutls.nakaRushtonDirections;
baseline              = p.Results.baseline;

% Create fit objects
% IAMP Object
if ~isempty(paramsIAMP) | ~isempty(paramsMeanBetas)
    temporalFitIAMPObj = tfeIAMP('verbosity','none');
end

% QCM Object
if ~isempty(paramsQCMDir)
    QCMDirectionObj = tfeQCMDirection('verbosity','none','dimension',analysisParams.theDimension);
end

% Naka-Ruston Object
if ~isempty(paramsNakaRushton)
    commonOffsetNRObj = tfeNakaRushtonDirection(nakaRushtonDirections);
end

% Get subplot sizing
rws = ceil(sqrt(length(packets)));
cols = rws-1;
if rws*cols < length(packets)
    cols = rws;
end

% Set indexing for betas
if ~isempty(paramsMeanBetas)
    betaLength = (length(paramsMeanBetas.paramMainMatrix)-1)/length(analysisParams.sessionFolderName);
end

% Open figure
figure

% Use QCM fit to IAMP to predict timecourse.
counter = 1;
for ii = 1:length(analysisParams.sessionFolderName)
    for jj = 1:analysisParams.numAcquisitions
        
        % Plot
        subplot(rws,cols,counter); hold on
        
        % Plot fMRI Time Course
        plot(packets{counter}.response.timebase,packets{counter}.response.values,'Color',[.5 0 0]);
        
        % Plot IAMP predictions to stimulus
        if exist('paramsIAMP')
            IAMPResponses = temporalFitIAMPObj.computeResponse(paramsFitIAMP{counter},packets{counter}.stimulus,packets{counter}.kernel);
            plot(IAMPResponses.timebase, IAMPResponses.values,'Color',[.1 .8 0]);
        end
        
        % Doctor up the parameters to use mean IAMP values and plot again
        if exist('paramsMeanBetas')
        IAMPResponsesMean = temporalFitIAMPObj.computeResponse(paramsFitIAMPMean,packets{counter}.stimulus,packets{counter}.kernel);
        plot(IAMPResponsesMean.timebase,IAMPResponsesMean.values,'Color',[0 0.1 .9]);
        end

        % Doctor up parameters to use the QCM fit to the mean IAMP
        paramsFitIAMPQCM = paramsFitIAMP{counter};
        paramsFitIAMPQCM.paramMainMatrix(1:end-1) = [fitResponseStructQCM.values(1+((ii-1)*betaLength):ii*betaLength), baseline]' ;
        IAMPResponsesQCM = temporalFitIAMPObj.computeResponse(paramsFitIAMPQCM,packets{counter}.stimulus,packets{counter}.kernel);
        plot(IAMPResponsesQCM.timebase,IAMPResponsesQCM.values,'Color',[0 0 0]);
        
        % Set axis labels
        ylabel('PSC')
        xlabel('Time (mS)')
        title(sprintf('Session %s, Run %s', num2str(ii), num2str(jj)))
        % Change line size
        set(findall(gca, 'Type', 'Line'),'LineWidth',1);
        counter = counter +1;
    end
end
legend('time course','IAMP fit',' Mean IAMP params', 'Mean QCM params')
end