function[maxContrastActual,contrastDirectionActual] = getSubjectActualContrast(validationType)
% Get the actual contrast values from median of the valiadtionfile in the
% experiment for the positive and negative arms of the modulation with
% either the 2 or 15 degree cone fundamentals. 
% 
% Inputs:
%   validationType          = String. Either '2DegPostive', '2DegNegative',
%                             '15DegPostive',or '15DegNegative'.
%
% Outputs:
%   maxContrastActual       = The actual max contrast per direction.
%   contrastDirectionActual = Vector components of a unit vector pointed
%                             in the stimulus direction. 

% KAS25
switch validationType
    case '2DegPostive'
        coneContrastL = [0.0819,0.4274,0.1373,-0.0040,0.0801,0.1858,0.1507,-0.0510];
        coneContrastM = [-0.0858,0.4341,-0.0011,0.2225,-0.0253,0.0821,0.3754,0.1209];
    case '2DegNegative'
        % Sign flipped
        coneContrastL = [0.0906,0.4312,0.1436,0.0007,0.0761,0.1885,0.1531,-0.0502];
        coneContrastM = [-0.0879,0.4342,0.0007,0.2207,-0.0350,0.0785,0.3712,0.1155];
    case '15DegPostive'
        coneContrastL = [0.0799,0.4152,0.1336,0.0038,0.0787,0.1794,0.1578,-0.0483];
        coneContrastM = [-0.0873,0.4143,-0.0032,0.2229,-0.0240,0.0793,0.3671,0.1211];
    case '15DegNegative'
        coneContrastL = [0.0864,0.4183,0.1382,0.0066,0.0754,0.1822,0.1588,-0.0485];
        coneContrastM = [-0.0896,0.4147,-0.0024,0.2196,-0.0346,0.0738,0.3648,0.1147];
end

[theta,maxContrastActual] =  cart2pol(coneContrastL,coneContrastM);


[unitContrastL,unitContrastM] = pol2cart(theta,1);
contrastDirectionActual = [unitContrastL;unitContrastM;zeros(size(unitContrastM))];

end
