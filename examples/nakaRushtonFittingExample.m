%% Naka Ruston fitting example

%% Generate the data:

% Naka Ruston Params
Rmax = 1;
sigma = 5;
n = 1;
offset = 0;

% contrast values
C  = [.2, .4,.6,.8,1];

% response values
R  = nakaRushton(C,sigma,n,Rmax, offset);

% Add noise 
noiseLevel = 0;
R_noise = R + rand(size(R)).*noiseLevel.*Rmax;

%% Fit the data 
%   Rmax = params(1)
%   sigma = params(2)
%   n = params(3)
%   response = Rmax*[contrast^n]/[contrast^n + sigma^n]
params0 = [];
[params,f] = FitNakaRushton(C,R_noise,params0)

% Get the naka predictions for the inital contrast values.
nrResponses = nakaRushton(C,params(2), params(3),params(1), offset);

% upsample the function
contrastSpacing = linspace(max(C),min(C),30);
nrResponses_upsampled = nakaRushton(contrastSpacing,params(2), params(3),params(1), offset);

%% Plot the results
figure; hold on
p1 = scatter(C,R_noise,'^k');
p2 = scatter(C,nrResponses,'*r');
p3 = plot(contrastSpacing,nrResponses_upsampled,'k--');
xlabel('Contrast')
ylabel('Response')
legend([p1 p2 p3], 'Data Points','NR Fits', 'NR Plot','Location','northwest')
