function response = nakaRushton(contrastSupport,K,n,Rmax,offset)
% Compute Naka-Rushton equation with offset
%
% Syntax:
%   response = nakaRushton(contrastSupport,K,n,Rmax,offset)
%
% Description: 
%   This function returns output of a contrast response function
%   as estimated by the Naka-Rushton function. Rmax controls the "percent
%   signal change", K is the the semisaturation constant, and n controls
%   the slope.
%  
%   response = Rmax*[contrastSupport^n]/[contrast^n + K^n] + offset
%
% Inputs: 
%   contrastSupport  - The contast value support for the naka rushton
%                      function. Either scalar or vector of contrast
%                      values.
%   K                - The semisaturation constant. 
%   n                - Controls the slope of the function. 
%   Rmax             - Maximun response value.
%   offset           - Function amplitude offset.
%
% Outputs: 
%   response         = Simulated response to the contrast values.
%
% See also: ComputeNakaRushton, InvertNakaRushton

% History:
%   03/08/17  mab    Wrote it.
%   11/24/18  dhb    Call through BLTB ComputeNakaRushton.

% Naka-Rushton Function
response = ComputeNakaRushton([Rmax,K,n],contrastSupport) + offset;
%response = Rmax .* ( contrastSupport.^n ./ (contrastSupport.^n + K.^n)) + offset;

end
