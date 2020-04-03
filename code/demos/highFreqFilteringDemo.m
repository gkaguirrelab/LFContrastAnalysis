%%Time specifications:
Fs = 1/.800;                 % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 288;              % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc = 1/StopTime;                     % hertz
x1 = cos(2*pi*Fc*t);

Fc = 2/StopTime;                     % hertz
x2 = cos(2*pi*Fc*t);

Fc = 3/StopTime;                     % hertz
x3 = cos(2*pi*Fc*t);

Fc = 24/StopTime;                     % hertz
x4 = cos(2*pi*Fc*t);

S = x1+x2+x3+x4 %+ randn(size(t))/2;

Y = fft(S);
L = length(Y);

figure;
subplot(2,2,1);
plot(t,S);
title('X(t)')
xlabel('Response')
ylabel('Time (s)')

subplot(2,2,2);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

sHat = highpass(S,12/288,1/.8,'Steepness',0.95); 
yHat = fft(sHat);
L = length(yHat);


subplot(2,2,3);
plot(t,sHat);
title('Highpass X(t)')
xlabel('Response')
ylabel('Time (s)')

subplot(2,2,4);
P2 = abs(yHat/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

plot(f,P1) 
title('Amplitude Spectrum of High Pass X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')




