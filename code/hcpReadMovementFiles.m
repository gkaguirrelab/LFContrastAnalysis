function [movenentRegressor] = hcpReadMovementFiles(txtFile,rws,cols)

fileID = fopen(txtFile,'r');
formatSpec = '%f';
textVector = fscanf(fileID,formatSpec);
fclose(fileID);

movementRegressorsFull = reshape(textVector,[cols,rws])';
movenentRegressor = movementRegressorsFull(:,7:12);

end

