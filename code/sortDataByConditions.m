function [avgPerCond, blockAvg] = sortDataByConditions(meanSignal,expParams)

for ii = 1:length(expParams)
    blockAvg(ii) = mean(meanSignal(expParams(ii,1):expParams(ii,2)));
end

conditions = unique(expParams(:,3));

for jj = 1:length(conditions)
    avgPerCond(jj) = mean(blockAvg(expParams(:,3) == conditions(jj))); 
end

end