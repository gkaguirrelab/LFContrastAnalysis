function LMSContrastMat = LMSContrastValuesFromParams(expParams,contrastCoding,directionCoding,maxContrastPerDir,totalTime,deltaT)
% LMSContrastValuesFromParams  Loop over trials, show stimuli and get responses.
%
% Syntax:
%    LMSContrastMat = LMSContrastValuesFromParams(expParams,contrastCoding,directionCoding,maxContrastPerDir,totalTime,deltaT)
%
% Description:
%    This function takes in a description of the stimulus order for a
%    given run and returns a 3xn matrix that contains the L, M and, S contrast
%    for each time point.
%
% Input:
%    expParams (mat)       - An nx4 martix that destribes the timing and
%                            order of the stimulus. The 1st column is the
%                            block start times, the 2nd column is the block
%                            end times, the 3 columns is the contast level
%                            code (this the code used in TrialSequenceMR
%                            to index the different contrast), the fourth
%                            column is color directrion code used in the
%                            experiment (ex. 1 = L+M, 2 = L-M,...).
%    contrastCoding (vec)  - Contrast levels that match the order of the
%                            contrast level code in expParams (ex. [100,
%                            50, 25, 12.5, 6.25, 0])
%    directionCoding (mat) - A matrix descibing the contrast coding for the
%                            color directions. matrix a 3xn matrix where n
%                            is the number of directions.
%    maxContrastPerDir(vec)- A vector in the same order as the direction coding
%                            that desribes the max contrast for that
%                            direction
%    totalTime (scalar)    - total time of run
%    deltaT (scalar)       - timebase resolution
%
% Output:
%    regressors
%
% Optional key/value pairs:
%    none

% Example:
% expParams = getExpParams(dataParamFile,TR,'hrfOffset', false, 'stripInitialTRs', false);
% contrastCoding = [1, .5, .25, .125, .0625, 0];
% directionCoding = [1,1,1,0;-1,1,0,1;0,0,0,0] %this 1 = L-M 2 = L+M 3 = L 4 = M; 
% maxContrastPerDir = [.6,.40,.10,.10] % max contrast in the same order as above
% totalTime = protocolParams.nTrials * protocolParams.trialDuration * 1000;
% deltaT = 800;
% LMSContrastValuesFromParams(expParams,contrastCoding,directionCoding,maxContrastPerDir,totalTime,deltaT)

LMSContrastMat = zeros(3,totalTime/deltaT);

for kk = 1:size(expParams,1)
    LMSbaseDir = repmat(directionCoding(:,expParams(kk,4)), [1,(expParams(kk,2)-expParams(kk,1)+1)]);
    LMSmaxContrast = LMSbaseDir.* maxContrastPerDir(expParams(kk,4));
    LMScontrastBlock = LMSmaxContrast .* contrastCoding(expParams(kk,3));
    LMSContrastMat(:,expParams(kk,1):expParams(kk,2)) = LMScontrastBlock;    
end

end
