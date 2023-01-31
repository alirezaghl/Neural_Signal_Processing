load EEGrestingState.mat
npnts = length(eegdata);
time = (0:npnts-1)/srate;

figure(1)
subplot(2,1,1)
plot(time,eegdata,'b')
xlabel('Time(s)'), ylabel('Voltage(\muv)')

%static specteral analysis
hz = linspace(0,srate/2,floor(npnts/2)+1);
amplitude = 2*abs(fft(eegdata)/npnts); 
power = amplitude.^2;
subplot(2,1,2)
plot(hz,amplitude(1:length(hz)),'k','linew',2)
hold on
plot(hz,power(1:length(hz)),'r','linew',2)
xlabel('Frequecny(hz)'), ylabel('Amplitude/Power')
legend('Amplitude','Power')

clc
clear
close all
load restingstate64chans.mat

% compute the power spectrum of each epoch seperately then average them
% together
eeg.data = double(EEG.data); % 64*2048*63 
chanpower = (2*abs(fft(eeg.data,[],2)/EEG.pnts)).^2;
chanpower = mean(chanpower,3);
hz = linspace(0,EEG.srate/2,floor(EEG.pnts/2)+1);
figure(2)
plot(hz,chanpower(:,1:length(hz)), 'linew',2)
set(gca,'xlim',[0 30], 'ylim', [0 50])

% Topographical map
alphabounds = [8 12];
freqidx = dsearchn(hz', alphabounds');
alphapower = mean(chanpower(:,freqidx(1):freqidx(2)),2);
topoplotIndie(alphapower,EEG.chanlocs)
colormap hot