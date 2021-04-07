function [lPlusMDist, lMinusMDist] = invertQCMforCardinal(targetResponse,qcmParams)
% This function returns the distane 

%% get the eq. contrast from NR params
eqContrast = InvertNakaRushton([qcmParams.crfAmp,qcmParams.crfSemi,qcmParams.crfExponent],targetResponse-qcmParams.crfOffset);

if isnan(eqContrast)
    lPlusMDist = nan;
    lMinusMDist = nan;
else
%% get the coeff of the ellipse
% scale matrix with the minor axis ratio param
Smat = [1, 0; 0, 1/qcmParams.Qvec(1)];
% rotation matrix with the ellipse angle param
Vmat = deg2rotm(qcmParams.Qvec(2));
% multiply to get the Q matrix which gives us the coeffs of the ellipse

Q = Vmat*Smat*Smat'*Vmat';
% splitting up the coeffs
a = Q(1,1);
b = Q(1,2);
c = Q(2,2);

%% Get the distance to the ellipse for the 45 and 315 direction
% The following equation is the polar equation for an ellipse: 
% plugging in r*cos(?) for x and r*sin(?) for y into: 
%    a*x^2 + 2*b*x*y + c*y^2 = ±eqContrast
% which gives:
%    a*r^2*cos^2(?) + 2*b*r^2*cos(?)*sin(?) + c*r^2*sin^2(?) = ±eqContrast
% Solve for r gets:
%    r = sqrt(eqContrast/(a*cosd(?)^2 + 2*b*cosd(?)*sind(?) + c*sind(?)^2))

theta = 45;
lPlusMDist = sqrt(eqContrast./(a*cosd(theta).^2 + 2*b*cosd(theta)*sind(theta) + c*sind(theta).^2)); 

theta = 315;
lMinusMDist = sqrt(eqContrast./(a*cosd(theta).^2 + 2*b*cosd(theta)*sind(theta) + c*sind(theta).^2)); 
end

end
