% Creating filter
srate = 1024;
nyquist = srate/2;
frange = [20 25];
transw = .1;
order = round(10*srate/frange(1));
shape = [0 0  1 1 0 0];
frex = [0 frange(1) - frange(1)*transw frange frange(2) + frange(2)*transw nyquist] /nyquist ;
filtkernel = firls(order,frex,shape);
filter_power = abs(fft(filtkernel)).^2;
hz = linspace(0,nyquist,floor(length(filter_power)/2)+1);
figure(1)
subplot(2,1,1)
plot(filtkernel,'k','LineWidth',2)
axis square
subplot(2,1,2)
plot(hz,filter_power(1:length(hz)),'ks-','LineWidth',2)
hold on
plot(frex*nyquist,shape,'ro-','LineWidth',2)
set(gca, 'xlim',[0 80])
xlabel('Frequency(Hz)'), ylabel('Filter Gain')
title('Frequency responses of filter')
axis square

% Apply filter to random noise
signal = randn(srate*4,1);
figure(2)
subplot(2,2,1)
plot(signal,'k',LineWidth=2)
set(gca,'xlim', [0 4000])
axis square
subplot(2,2,2)
hz = linspace(0,nyquist,floor(length(signal)/2)+1);
sig_woutfilter = abs(fft(signal)).^2;
plot(hz,sig_woutfilter(1:length(hz)),'k',LineWidth=2)
axis square
set(gca,'xlim', [0 40])
subplot(2,2,3)
filtsignal = filtfilt(filtkernel,1,signal);
plot(filtsignal,'r',LineWidth=2)
set(gca,'xlim', [0 4000])
axis square
sig_wfilter = abs(fft(filtsignal)).^2;
subplot(2,2,4)
plot(hz,sig_wfilter(1:length(hz)),'r','LineWidth',2)
set(gca,'xlim', [0 40])
axis square

% Hilber transform 
Hilbfiltsignal = hilbert(filtsignal);
figure(3)
subplot(3,1,1)
plot(real(Hilbfiltsignal))
% Magnitude
subplot(3,1,2)
plot(abs(Hilbfiltsignal))
% Phase
subplot(3,1,3)
plot(angle(Hilbfiltsignal))