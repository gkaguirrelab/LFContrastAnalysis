function [angles, quadrant] = getAnglesFromContrastMeasurements(contrastMat)


L(1) = contrastMat(1,1);
L(2) = contrastMat(1,2);
M(1) = contrastMat(2,1);
M(2) = contrastMat(2,2);

% Calculate positive and negative modulation angle (M = y & L = x) 
for ii = 1:length(L)
    
    % No Angle
    if M(ii) == 0 && L(ii) == 0
        angles(ii) = NaN;
        quadrant(ii)   = NaN;
    % Quadrant I
    elseif M(ii) >=0  && L(ii) >= 0
        angles(ii) = atand(abs(M(ii)/L(ii)));
        quadrant(ii) = 1;
    % Quadrant II 
    elseif M(ii) >= 0 && L(ii) <= 0
        angles(ii) = atand(abs(L(ii)/M(ii))) + 90;
        quadrant(ii) = 2;  
    % Quadrant III
    elseif M(ii) <= 0 && L(ii) <= 0
        angles(ii) = atand(abs(M(ii)/L(ii))) + 180;
        quadrant(ii) = 3;
    % Quadrant IV
    elseif M(ii) <= 0 && L(ii) >=0
        angles(ii) = atand(abs(L(ii)/M(ii))) + 270; 
        quadrant(ii) = 4;
    % Weird Stuff
    else
        angles(ii) = NaN;
        quadrant(ii) = NaN;
    end
        
end