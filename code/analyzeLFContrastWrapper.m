% THis is a script to loop over subjects and generate figures

% Cell array fo subject names ,'KAS25','KAS25_replication'
subjIds = {'LZ23','LZ23_replication','AP26','AP26_replication'};

% loop over subjects

for ii = 1:length(subjIds)
    % Make the CRF, time couse and, the ellipse non-lin plots
    analyzeLFContrast(subjIds{ii})
    
    % Make the cross validation time course and R^2 plots
    analyzeLFContrast_crossValidate(subjIds{ii})
end

subjIds = {'KAS25','KAS25_replication','LZ23','LZ23_replication','AP26','AP26_replication'};
for jj = 1:length(subjIds)
    % Make Maps of parameters
    analyzeLFContrast_voxelwise(subjIds{jj})
end


% Make scatter plots
subjIds = {'KAS25','KAS25_replication','LZ23','LZ23_replication','AP26','AP26_replication'};
maps = {'minorAxis','angle','rSqaured'};
for kk = 1:length(subjIds)
    for ll = 1:length(maps)
        scatterMapWithEcc(subjIds{kk}, maps{ll});
    end
end