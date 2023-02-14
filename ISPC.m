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