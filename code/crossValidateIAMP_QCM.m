function [rmseMeanIAMP rmseMeanQCM rmseSemIAMP rmseSemQCM] = crossValidateIAMP_QCM(betas)
% Calculates the cross validated RMSE and SEM for both the QCM and IAMP  
%
% Syntax:
%
%
% Description:
%  
%
% Inputs:
%
%
% Outputs:
%    rmseMeanIAMP
%    rmseMeanQCM
%    rmseSemIAMP
%    rmseSemQCM
%
% Optional key/value pairs:
%    none

% Pick runs to leave out
leaveOut =[];
for ii = 1:size(betas,3)
    leaveOut = [leaveOut; randperm(size(betas,2))];
end



for jj = 1:size(betas,2)
    for kk = 1:size(leaveOut,1)
        sheet = betas(:,:,kk);
        sheet(:,leaveOut(kk,jj)) =[];
        betasCV(:,:,kk) = sheet;
    end
    
    baselineCV = betasCV(end,:,:);
    betasCV = betasCV(1:end-1,:,:);
    
    tmpBetas = [];
   for kk = 1:size(betasCV,3)
      tmpBetas =[tmpBetas; betasCV(:,:,kk)]; 
   end
   
   meanBetasCV = [mean(tmpBetas,2); mean(baselineCV(:))];
       
   
   clear betasCV baselineCV
end

