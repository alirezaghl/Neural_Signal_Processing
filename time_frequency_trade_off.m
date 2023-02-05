srate = 512;
num_cycles = [2 6 8 15];
frex = 8;
time = -2:1/srate:2;
hz = linspace(0,srate/2,floor(length(time)/2)+1);

for i=1:4
    subplot(4,2,i*2-1)
    s = num_cycles(i) / (2*pi*frex);
    plot(time,exp((-time.^2)./(2*s^2) ),'k','LineWidth',3)
    title([ 'Wavelet at ' num2str(num_cycles(i)) ' cycles' ])
    % frequency domain
    subplot(4,2,i*2)
    cmw = exp(1i*2*pi*frex.*time) .* exp((-time.^2)./(2*s.^2));
    cmwX = fft(cmw);
    cmwX = cmwX./max(cmwX);
    plot(hz,abs(cmwX(1 : length(hz))),'k','LineWidth',3)
    set(gca,'xlim',[0 20])
end

clear
load sampleEEGdata.mat
num_frex = 40;
min_frequency = 2;
max_frequency = 40;
channel2use = 'o1';
num_cycles = [2 6 8 15];
baseline_window = [-500 -200];
frex = linspace(min_frequency,max_frequency,num_frex);
time = -2:1/EEG.srate:2;
half_wave = (length(time)-1)/2;
nkern = length(time);
ndata = EEG.pnts * EEG.trials;
nconv = nkern  + ndata -1;
tf = zeros(length(num_cycles),length(frex),EEG.pnts);
baseidx = dsearchn(EEG.pnts',baseline_window');
eegfft = fft(reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),nconv);

for cyclei=1:length(num_cycles)
    for fi=1:length(frex)
        s = num_cycles(cyclei) / (2*pi*frex(fi));
        cmw = exp(1i * 2 * pi * frex(fi) .* time) .* exp((-time.^2)./(2*s^2));
        cmwX = fft(cmw,nconv);
        cmwX = cmwX ./max(cmwX);
        as = ifft(cmwX.*eegfft);
        as = as(half_wave+1:end-half_wave);
        as = reshape(as,EEG.pnts,EEG.trials);
        tf(cyclei,fi,:) = mean(abs(as).^2,2);
    end
    tf(cyclei,:,:) = 10*log10(bsxfun(@rdivide,squeeze(tf(cyclei,:,:)),mean(tf(cyclei,:,baseidx(1):baseidx(2)),3)' ) ); %baseline normalization
end

figure

for cyclei=1:length(num_cycles)
    subplot(2,2,cyclei)
    contourf(EEG.times,frex,squeeze(tf(cyclei,:,:)),40,'linecolor','none')
    set(gca,'clim',[-3 0] * 2, 'ydir', 'normal','xlim',[-300 1000])

end

clc
clear
close all

range_cycles = [4 13];
ncycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex); % creating range of cycles 
tf = zeros(length(frex),EEG.pnts);
for fi=1:length(frex)
        s = ncycles(fi) / (2*pi*frex(fi));
        cmw = exp(1i * 2 * pi * frex(fi) .* time) .* exp((-time.^2)./(2*s^2));
        cmwX = fft(cmw,nconv);
        cmwX = cmwX ./max(cmwX);
        as = ifft(cmwX.*eegfft);
        as = as(half_wave+1:end-half_wave);
        as = reshape(as,EEG.pnts,EEG.trials);
        tf(fi,:) = mean(abs(as).^2,2);
end

tfDB = 10*log10(bsxfun(@rdivide,tf,mean(tf(:,baseidx(1):baseidx(2)),2) ) ); % baseline normalization
figure(2)
subplot(2,1,1)
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 5] , 'ydir', 'normal','xlim',[-300 1000])
colormap hot
colorbar

subplot(2,1,2)
contourf(EEG.times,frex,tfDB,40,'linecolor','none')
set(gca,'clim',[-3 3] , 'ydir', 'normal','xlim',[-300 1000])
colormap hot
colorbar