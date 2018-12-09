%% This demo is designed to simulte the LFContrast analysis pipeline

%--------------------------------------------------------------------------
% Step 1:  Set up general params.
%--------------------------------------------------------------------------

% the directions in the LM plane of the stimulus vector
angles =[-45,0,45,90];
contrastDirections = vectorAngle2LMScontrast(angles,'LM');
numDirections      = size(contrastDirections,2);

% Contrast values presented per direction.
contrastLevels     = [1.0, 0.5, 0.25, 0.125, 0.0625];
numContrastLevels  = length(contrastLevels);

% Add noise
noiseLevel = .2;

% QCM ellipse number of points 
nQCMPoints = 200;

%--------------------------------------------------------------------------
%% Step 2:  Generate the data.
%--------------------------------------------------------------------------
for ii = 1:numDirections
    % Naka Ruston Params
    Rmax   = (1.35-.60).*rand(1) + .60;
    sigma  = (.75-.35).*rand(1) + .35;
    n      = 2;
    offset = 0;
    
    % response values
    R  = tfeQCMComputeNakaRushton(contrastLevels,sigma,n,Rmax, offset);
    betas(ii,:) = R + rand(size(R)).*noiseLevel.*Rmax;
end
tmp = betas';
allBetas = tmp(:);

% Plot the data
figure; hold on
for jj = 1:numDirections
    plot(contrastLevels,betas(jj,:))
end
% concat betas
betas

%--------------------------------------------------------------------------
%% Step 3:  Fit the 
%--------------------------------------------------------------------------

% Fit the "betas" with QCM
% Set parameters and construct a QCM object.
temporalFitQCM = tfeQCM('verbosity','none','dimension',2);

% Set up contrast values matched to resoponse order
% Set up stim order info to creat LMS contrast by timepoint matrix
stimulusStruct.values   = generateStimCombinations(contrastLevels,contrastDirections,[1,1,1,1],2);
stimulusStruct.timebase = 1:length(stimulusStruct.values);

% Snag response values from IAMP fit.
%end -1 is for the attention event modeling
thePacket.response.values = allBetas';
thePacket.response.timebase = 1:length(thePacket.response.values);

% Construct a packet for the QCM to fit.
thePacket.stimulus = stimulusStruct;
thePacket.kernel = [];
thePacket.metaData = [];

% Fit
[paramsQCMFit,fVal,fitResponseStructQCM] = temporalFitQCM.fitResponse(thePacket);
fprintf('Model parameter from fits:\n');
temporalFitQCM.paramPrint(paramsQCMFit)





eqContrast = InverttfeQCMComputeNakaRushton([paramsQCMFit.crfAmp,paramsQCMFit.crfSemi,paramsQCMFit.crfExponent],thresh);
circlePoints = eqContrast*UnitCircleGenerate(nQCMPoints);
[~,Ainv,Q] = EllipsoidMatricesGenerate([1 paramsQCM.Qvec],'dimension',2);
ellipsePoints = Ainv*circlePoints;
checkThresh = ComputetfeQCMComputeNakaRushton([paramsQCM.crfAmp,paramsQCM.crfSemi,paramsQCM.crfExponent],diag(sqrt(ellipsePoints'*Q*ellipsePoints)));
if (any(abs(checkThresh-thresh) > 1e-10))
    error('Did not invert QCM model correctly');
end

%% Plot data points
if (isempty(hdl))
    hdl = figure; hold on
else
    figure(hdl); hold on
end
sz = 50;
scatterHdl = scatter(dataPointsNR(:,1),dataPointsNR(:,2),sz,'MarkerEdgeColor',color,'MarkerFaceColor',color,'LineWidth',1.5);
scatter(dataPointsLI(:,1),dataPointsLI(:,2),sz,color,'x')


% Add ellipse
plot(ellipsePoints(1,:),ellipsePoints(2,:),'color', color);





