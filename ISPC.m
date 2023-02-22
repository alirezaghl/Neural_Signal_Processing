load v1_laminar.mat
chanl1idx = 1;
chanl2idx = 8;
freq = 8;
% creating complex morlet wavelet
time = -2:1/srate:2;
s = 8/(2*pi*freq);
cmw = exp(1i * 2 * pi * freq .* time) .* exp(-time.^2./(2*s^2));
figure(1)
plot(time,cmw)
half_wave = (length(time)-1)/2;
nkern = length(time);
ndata = size(csd,2);
nconv = ndata + nkern - 1;
% fourier transform of wavelet
waveletx = fft(cmw,nconv);
waveletx = waveletx./max(waveletx);
figure(2)
plot(abs(waveletx))
% fourier transform of channel_1 data
phase_data = zeros(2,length(timevec));
real_data = zeros(2,length(timevec));
ch1x = fft(csd(1,:,1),nconv);
ch1_ifft = ifft(ch1x.*waveletx,nconv);
ch1_ifft = ch1_ifft(half_wave+1:end-half_wave);
phase_data(1,:) = angle(ch1_ifft);
real_data(1,:) = real(ch1_ifft);
% fourier transform of channel_8 data
ch8x = fft(csd(8,:,1),nconv);
ch8_ifft = ifft(ch8x.*waveletx,nconv);
ch8_ifft = ch8_ifft(half_wave+1:end-half_wave);
phase_data(2,:) = angle(ch8_ifft);
real_data(2,:) = real(ch8_ifft);
% plot phase and real part of two channels
figure(3)
plot(timevec,phase_data)
figure(4)
plot(timevec,real_data)
% Phase synchronization
phase_synchronization = abs(mean(exp(1i*diff(phase_data))));

%% Phase lag index and ISPC on EEG data

load sampleEEGdata.mat
chan1 = 'FCz';
chan2 = 'POz';

min_freq = 2;
max_freq = 40;
num_frex = 34;
fwhm = linspace(.3,.1,num_frex);
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -2:1/EEG.srate:2-1/EEG.srate;
half_wave = (length(time)-1)/2;
nkern = length(time);
ndata = EEG.pnts * EEG.trials;
nconv = nkern  + ndata -1;
fft_ch1 = fft(reshape(EEG.data(strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),nconv);
fft_ch2 = fft(reshape(EEG.data(strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),nconv);
ispc = deal(zeros(num_frex,EEG.pnts));
PLI = deal(zeros(num_frex,EEG.pnts));

for fi=1:length(frex)
    cmw = exp(1i * 2 * pi * frex(fi) .* time) .* exp((-4 * log(2) * time.^2)./(fwhm(fi)^2));
    cmwx = fft(cmw,nconv);
    cmwx = cmwx ./ max(cmwx);
    as1 = ifft(cmwx.*fft_ch1);
    as1 = as1(half_wave+1:end-half_wave);
    as1 = reshape(as1,EEG.pnts,EEG.trials);
    as2 = ifft(cmwx.*fft_ch2);
    as2 = as2(half_wave+1:end-half_wave);
    as2 = reshape(as2,EEG.pnts,EEG.trials);
    cdd = exp(1i*( angle(as1)-angle(as2)));
    ispc(fi,:)  = abs(mean(cdd,2));
    PLI(fi,:)   = abs(mean(sign(imag(cdd)),2));
end

figure(1)
subplot(221)
contourf(EEG.times,frex,ispc,40,'linecolor','none')
set(gca,'xlim',[-300 1200],'clim',[0 .5])
colormap hot; colorbar
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('ISPC, voltage') 

subplot(222)
contourf(EEG.times,frex, PLI, 40, 'linecolor','none')
set(gca,'xlim',[-300 1200],'clim',[0 .5])
colormap hot; colorbar
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('PLI, voltage')
