function contrast = vectorAngle2LMScontrast(angles, plane, varargin)
% This function returns the vector components of the contrast within a 2
% dimensional plane in cone space
%
% Syntax:
%   contrast = vectorAngle2LMScontrast(angles, plane, varargin)
%
% Description:
%    This function plots the IAMP fits and IAMP-QCM predictions from 
%    runIAMP_QCM.m as contrast response functions (one plot per modulation 
%    direction)
%
% Inputs:
%    angles       - input angle(s) in the 2 dimensional cone contrast plane 
%                   (scalar or vector)
%    
%    plane        - Set the X and Y axis of the space as L,M, or S. (only 
%                   the LM is coded as of now.)
%
% Outputs:
%    contrast     - Cone contrast components of the input contrast in the 
%                   desired space 
%
% Optional key/value pairs:
%    precision    - Precision for rounding the output contrasts (default= 4)

% MAB 09/18

p = inputParser;
p.addRequired('angles',@isnumeric);
p.addRequired('plane',@ischar);
p.addParameter('precision',4,@isnumeric);
p.parse(angles, plane, varargin{:});

for ii = 1:length(angles)
    switch plane
        case 'LM'
            contrast(1,ii) = cosd(angles(ii));
            contrast(2,ii) = sind(angles(ii));
            contrast(3,ii) = 0;
        case 'MS'
            % this needs to be though about more. only need LM plane now.
        case 'LS'
            % this needs to be though about more. only need LM plane now.
    end
end

contrast = round(contrast,p.Results.precision);
end
