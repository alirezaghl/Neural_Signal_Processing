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





