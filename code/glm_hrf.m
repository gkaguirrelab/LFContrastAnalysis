function hrf = glm_hrf(t, varargin)
%
%    hrf = boyntonHIRF(t, [n=3], [tau=1.08], [delay=2.05]])
%
%Purpose:
%   Compute the Boynton et al. HIRF function.  The return values include
% both the hrf and the values of time and the parameters used in the
% computation.
%
% Equation:  (Eq 3 from Boynton & Heeger, J Neurosci 1996)
%   h(t) = [(t/tau) ^ (n-1) * exp(-t/tau)] / [tau(n-1)!]
%
% Inputs:
%   t: t should be in units of SECONDS, and reflect the time window in
%   which to estimate the HRF. 
%
%   n: exponent in eq. above. (corresponds to z in standard gamma
%   functions).
%
%   tau: time constant in eq.
%
%   delay: additional delay before onset of gamma function. This is an
%   added heuristic, which may vary from subject to subject.
%
%
% Outputs:
%   hrf: estimate of the HRF sampled at t seconds.
%   IMPORTANT: If you are going to convolve
%   this for fMRI analysis, and the TR of your data is not 1, you will
%   need to resample both the input t and hrf to match the MR frames.
%
% written 2005 by wandell.
% ras, 01/2007: heavily modified: no parms struct, each arg can be
% specified separately; clarifies difference between seconds and frames,
% doesn't modify t.

%Input Parser
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;
p.addRequired('t',@isvector);
p.addParameter('n',3,@isnumeric);
p.addParameter('tau',1.08,@isnumeric)
p.addParameter('delay',2.05,@isnumeric);
p.parse(t,varargin{:});

% Set up params for main computation 
n     = p.Results.n;
tau   = p.Results.tau;
delay = p.Results.delay;

% initialize the HRF to be zeros, the same size as t:
hrf = zeros( size(t) );

% The HRF is not specified for t < 0 secs. In addition, any
% values before the delay should also be zero.
% So, we only sample the HRF below, for time points after [0+delay].
% We call this sampling vector x to distinguish it from t.
x = t - delay;
x = x(x>0);  % not defined for x<0

if round(n)~=n
    warning('Boynton HIRF parameters stored as unusual values.  n is being rounded from %.2f to %i.',n, round(n))
    disp('Check stored glmHRF_params against cannonical values.')
    n = round(n);    
end

% main computation (per Boynton & Heeger, 1996)
tmpHrf = (x/tau).^(n-1) .* exp(-(x/tau)) / (tau*(factorial(n-1)));

% paste in tmpHRF into appropriate indices in hrf, corresponding
% to the (shifted after delay) time points:
hrf(end-length(tmpHrf)+1:end) = tmpHrf;

return;