function regressors = createRegressors(expParams,baselineCondNum,totalTime,deltaT)


condRegMat = zeros(length(unique(expParams(:,3)))-length(baselineCondNum),totalTime/deltaT,length(unique(expParams(:,4))));
baselineRegVec = zeros(1,totalTime/deltaT);

for kk = 1:size(expParams,1)
    if expParams(kk,3) ~= baselineCondNum
        condRegMat(expParams(kk,3),expParams(kk,1):expParams(kk,2),expParams(kk,4)) = 1;
    elseif expParams(kk,3) == baselineCondNum
        baselineRegVec(1,expParams(kk,1):expParams(kk,2)) = 1;
    end
end
stimulusStruct.values = [];
for ii = 1:length(unique(expParams(:,4)))
    stimulusStruct.values = vertcat(stimulusStruct.values,condRegMat(:,:,ii));
end

regressors = vertcat(stimulusStruct.values,baselineRegVec);

end