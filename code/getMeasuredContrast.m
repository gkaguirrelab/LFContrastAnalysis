for jj = 1:4
    size(directedDirection{jj}.describe.validation)
    count = 1;
    for ii = 6:15
        contrasts(:,count) = mean(abs(directedDirection{jj}.describe.validation(ii).contrastActual),2);
        count = count +1;
    end
    
    contrastsCenter = contrasts(1:3,:);
    contrastsPeriph = contrasts(4:6,:);
    
    contCent = vecnorm(contrastsCenter(1:2,:));
    contPeriph = vecnorm(contrastsPeriph(1:2,:));
    
    meanContCent(jj) = mean(contCent);
    meanContPeriph(jj) = mean(contPeriph);
    
    stdContCent(jj) = std(contCent,0,2);
    stdContPeriph(jj) = std(contPeriph,0,2);
    
    anglesCent = atand(contrastsCenter(2,:)./contrastsCenter(1,:));
    anglesPeriph =  atand(contrastsPeriph(2,:)./contrastsPeriph(1,:));
    
    meanAnglesCent(jj) = mean(anglesCent);
    meanAnglesPeriph(jj) = mean(anglesPeriph);
    
    stdAnglesCent(jj) = std(anglesCent);
    stdAnglesPeriph(jj) = std(anglesPeriph); 
end