%% Time frequency analysis of v1 dataset by wavelets
% the first problem is to extract power and phase from trial 10 and plot
% them
load v1_laminar.mat
num_frex = 43;
min_frequency = 10;
max_frequency = 100;
frex = linspace(min_frequency,max_frequency,num_frex);
time = -2:1/srate:2-1/srate;
half_wave = (length(time)-1)/2;
nkern = length(time);
ndata = length(timevec);
nconv = ndata+nkern-1;
range_cycles = [3 15];
ncycles = linspace(range_cycles(1),range_cycles(end),num_frex);
cmwx = zeros(length(frex),nconv); % size = 43 * 4577

for fi=1:length(frex)
        s = ncycles(fi) / (2*pi*frex(fi));
        cmw = exp(1i * 2 * pi * frex(fi) .* time) .* exp((-time.^2)./(2*s^2));
        cmwx(fi,:) = fft(cmw,nconv);
        cmwx(fi,:) = cmwx(fi,:)./max(cmwx(fi,:));
end

tf = zeros(size(csd,1), length(frex),length(timevec),2);

for chani=1:size(csd,1)
    lfpx = fft(squeeze(csd(chani,:,10)),nconv); % size = 1 * 4577
    for fi=1:length(frex)
        as = ifft(cmwx(fi,:).*lfpx);
        as = as(half_wave+1:end-half_wave);
        as_power = abs(as).^2;
        as_phase = angle(as);
        tf(chani,fi,:,1) = as_power;
        tf(chani,fi,:,2) = as_phase;

    end

end
figure(1)
subplot(2,1,1)
contourf(timevec,frex,squeeze(tf(6,:,:,1)),40,'linecolor','none')
set(gca, 'xlim', [-.2 1] ,'clim',[0 90000])
xlabel('Time(s)'), ylabel('Frequencies')
title('Power plot at channel 6')
axis square
subplot(2,1,2)
contourf(timevec,frex,squeeze(tf(6,:,:,2)),40,'linecolor','none')
set(gca, 'xlim', [-.2 1] ,'clim',[-3 3])
xlabel('Time(s)'), ylabel('Frequencies')
title('Phase plot at channel 6')
axis square

% the second problem is to extract power and phase from all trials and plot
% them
clear
load v1_laminar.mat
num_frex = 43;
min_frequency = 10;
max_frequency = 100;
frex = linspace(min_frequency,max_frequency,num_frex);
signal = csd(:,:,10);
time = -2:1/srate:2-1/srate;
half_wave = (length(time)-1)/2;
nkern = length(time);
ndata = length(timevec)*size(csd,3);
nconv = ndata+nkern-1;
range_cycles = [3 15];
ncycles = linspace(range_cycles(1),range_cycles(end),num_frex);
cmwx = zeros(length(frex),nconv);
for fi=1:length(frex)
        s = ncycles(fi) / (2*pi*frex(fi));
        cmw = exp(1i * 2 * pi * frex(fi) .* time) .* exp((-time.^2)./(2*s^2));
        cmwx(fi,:) = fft(cmw,nconv);
        cmwx(fi,:) = cmwx(fi,:)./max(cmwx(fi,:)); % size = 43 * 30845
end
tf = zeros(size(csd,1), length(frex),length(timevec),2);
lfpfft = fft(reshape(csd,size(csd,1), size(csd,2)*size(csd,3)),nconv,2); % 16 * 30845, 2 at the end -~ to take fft from timepnts
                                                                                        
% lfpx = fft(squeeze(csd(chani,:,:)),nconv);

for fi=1:length(frex)
        as = ifft(repmat(cmwx(fi,:),size(csd,1),1).*lfpfft,nconv,2);
        as = as(:,half_wave+1:end-half_wave);
        as = reshape(as,size(csd));
        as_power = mean(abs(as).^2,3);
        as_phase = abs(mean(exp(1i*angle(as)),3));
        tf(:,fi,:,1) = as_power;
        tf(:,fi,:,2) = as_phase;
end
figure(2)
subplot(2,1,1)
contourf(timevec,frex,squeeze(tf(6,:,:,1)),40,'linecolor','none')
set(gca, 'xlim', [-.2 1] ,'clim',[0 60000])
xlabel('Time(s)'), ylabel('Frequencies')
title('Power plot at channel 6')
axis square
subplot(2,1,2)
contourf(timevec,frex,squeeze(tf(6,:,:,2)),40,'linecolor','none')
set(gca, 'xlim', [-.2 1] ,'clim',[0 .25])
xlabel('Time(s)'), ylabel('Frequencies')
title('Phase plot at channel 6')
axis square