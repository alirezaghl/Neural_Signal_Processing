freq = 2;
srate = 1000;
time = -1:1/srate:1;
ampl = 2;
phas = pi/3;

sine_wave = ampl*sin(2*pi*freq*time+phas);
figure(1);
plot(time,sine_wave);

%% create multiple sine waves 
frequencies = [3 10 5 15 35];
amplit = [5 15 10 5 7];
phases = [pi/7 pi/8 pi pi/2 -pi/4];

sine_waves = zeros(length(frequencies), length(time)); % there are five same size sine waves

for fi=1:length(frequencies)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*frequencies(fi)*time+phases(fi));
end

% Plot sum of these five sine waves
figure(2), clf
plot(time, sum(sine_waves));

% Plot them seperately
figure(3), clf

for fi=1:length(frequencies)
    subplot(length(frequencies),1,fi);
    plot(time,sine_waves(fi,:))
end
