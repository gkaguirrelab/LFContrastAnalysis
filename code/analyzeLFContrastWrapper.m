% THis is a script to loop over subjects and generate figures 

% Cell array fo subject names
subjIds = {'KAS25','KAS25_replication','LZ23','LZ23_replication','AP26','AP26_replication'};

% loop over subjects
for ii = 1:length(subjIds)
    
    % Make the CRF, time couse and, the ellipse non-lin plots
    analyzeLFContrast(subjIds{ii})
    
    % Make the cross validation time course and R^2 plots
    analyzeLFContrast_crossValidate(subjIds{ii})
    
    % Make Maps of parameters
    analyzeLFContrast_voxelwise(subjIds{ii})
    
end
