% Creating a wavelet from scratch!
srate = 1000;
time = -3:1/srate:3;
frequecny = 2*pi;
ampl = 2;
phas = pi/2;

Cosine_wave = ampl*sin(2*pi*frequecny *time+phas);
figure
subplot(4,1,1)
plot(time,Cosine_wave,'k')

fwhm = .25; % width of the gausssian in seconds
gaussian_win = exp(-time.^2./(fwhm^2));
subplot(4,1,2)
plot(time,gaussian_win,'k')

morletW = gaussian_win.*Cosine_wave;
subplot(4,1,3)
plot(time,morletW,'k')

% Frequency domain
pnts = length(morletW);
morletW_amplitude = 2 * abs(fft(morletW,pnts)/pnts);
hz = linspace(0,srate/2,floor(pnts/2)+1);
morletW_power = morletW_amplitude.^2;
subplot(4,1,4)
plot(hz,morletW_power(1:length(hz)),'k')
set(gca,"XLim",[0,20])
xlabel('Frequecny(hz)'), ylabel('Power')

clc
clear
close all

% Complex Wavelet
sine_wave = exp(1i * 2 * pi * frequecny * time);
fwhm = .25; 
gaussian_win = exp(-time.^2./(fwhm^2));
complex_mw = gaussian_win .* sine_wave;
pnts = length(sine_wave);
morletW_amplitude = 2 * abs(fft(complex_mw,pnts)/pnts);
hz = linspace(0,srate/2,floor(pnts/2)+1);
morletW_power = morletW_amplitude.^2;
% plot(hz,morletW_power(1:length(hz)),'k','linew',2)
plot3(time,real(complex_mw),imag(complex_mw),'k','LineWidth',3)
rotate3d on

clc
clear 
close all
% time-frequency analysis of v1 data
figure(1)
load v1_laminar.mat
data = csd(6,:,10);
time = -1:1/srate:1;
subplot(6,1,1)
plot(timevec,data,'k')
set(gca,"XLim",[-.5,1.5])
frex = 45;
s=4/(2*pi*frex);
cmw = exp(1i * 2 * pi * frex * time) .* exp(-time.^2./(2*s^2)); % complex wavelet

nsignal = length(data);
nkernel = length(cmw);
nconvolution = nsignal + nkernel -1;
trim = floor(nkernel/2);
datax = fft(data,nconvolution);
hz = linspace(0,srate/2,floor(nconvolution/2)+1);
cmwx = fft(cmw,nconvolution);
subplot(6,1,2)
plot(hz,abs(datax(1:length(hz))),'k','LineWidth',2)
subplot(6,1,3)
plot(hz,abs(cmwx(1:length(hz))),'k','LineWidth',2)
cmwx = cmwx ./ max(cmwx); % Normalizing 
convresx = datax .* cmwx; % Creating the convolved signal with element-wise multiplication
subplot(6,1,4)
plot(hz,abs(convresx(1:length(hz))),'k'); % Ploting amplitude
convres = ifft(convresx);
convres = convres(trim+1:end-trim+1);
subplot(6,1,5)
plot(timevec,convres,'b')
set(gca,"XLim",[-.5,1.5])
subplot(6,1,6)
plot(timevec,data,'b')
hold on
plot(timevec,convres,'r')
legend('LFP data', 'Convolved data')
set(gca,"XLim",[-.5,1.5])

figure (2)
% plot imaginary part 
subplot(3,1,1)
plot(timevec, imag(convres))
% plot power
subplot(3,1,2)
plot(timevec, abs(convres).^2)
%plot phase
subplot(3,1,3)
plot(timevec, angle(convres))

% time-frequency analysis of v1 data, all trials!
% I used the variables from the last section
clc
clear
close all
load v1_laminar.mat

data = squeeze(csd(6,:,:));
dataR = reshape(data,1,[]); % reshape data from 1527 * 200 to 1 * 305400

nsignal = length(dataR);
nkernel = length(time);
nconvolution = nsignal + nkernel -1;
trim = floor(nkernel/2);
datax = fft(dataR,nconvolution);
cmwx = fft(cmw,nconvolution);
hz = linspace(0,srate/2,floor(nconvolution/2)+1);
convresx = datax .* cmwx;
convres = ifft(convresx);
convres = convres(trim+1:end-trim+1);
convresR = reshape(convres,1527,200);

figure(1)
subplot(1,2,1)
imagesc(timevec,[],data)
xlabel('Time'), ylabel('Trials')
set(gca, 'clim',[-1 1] * 2000)
subplot(1,2,2)
imagesc(timevec,[],abs(convresR))
xlabel('Time'), ylabel('Trials')
set(gca, 'clim',[-1 1] * 2000)
% plot power
figure(2)
plot(timevec, mean(abs(convresR).^2,2))

% time-frequency analysis of v1 data, all trials, with 30 diffrent
% frequency!
clc
clear
close all
load v1_laminar.mat

min_freq = 5;
max_freq = 90;
num_freq = 30;
frex = linspace(min_freq,max_freq,num_freq);
tf = zeros(num_freq,length(timevec));

for f=1:num_freq
    cmw = exp(1i * 2 * pi * frex(f) * time) .* exp(-4*log(2)*time.^2/.4^2); 
    cmwX = fft(cmw,nconvolution);
    cmwX = cmwX ./max(cmwX);
    as = ifft(datax .* cmwX);
    as = as(trim+1:end-trim+1);
    as = reshape(as,size(data));
    aspow = abs(as).^2;
    tf(f,:) = mean(aspow,2);
end    
contourf(timevec,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 1] * 1000, 'xlim', [-.1 1.4])
figure(2)
plot(timevec,mean(aspow,2))