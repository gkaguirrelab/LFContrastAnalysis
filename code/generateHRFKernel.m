function kernel = generateHRFKernel(gamma1,gamma2,gammaScale,timebase)
% Generates an HRF kernel for the TFE. 
%
% Syntax:
%   kernel = generateHRFKernel(gamma1,gamma2,gammaScale,timebase)
%
% Description:
%    This function takes in gamma functions parameters and returns a
%    difference of gammas HRF kernel. 
%
% Inputs:
%    gamma1            - positive gamma parameter (roughly, time-to-peak in
%                        secs).
%    gamma2            - negative gamma parameter (roughly, time-to-peak in
%                        secs).
%    gammaScale        - caling factor between the positive and negative
%                        gamma components.
%    timebase          - timebase for the kernel.
%
% Outputs:
%    kernel            - HRF kernel that can be used with the TFE.
%
% Optional key/value pairs:
%    none

% MAB 09/12/18

% The timebase is converted to seconds within the function, as the gamma
% parameters are defined in seconds.
hrf = gampdf(timebase/1000, gamma1, 1) - gampdf(timebase/1000, gamma2, 1)/gammaScale;

% prepare this kernelStruct for use in convolution as a BOLD HRF
kernelStruct.values   = hrf-hrf(1);
kernelStruct.timebase = timebase;

kernel = normalizeKernelArea(kernelStruct);
end