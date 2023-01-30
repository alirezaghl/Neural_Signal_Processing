srate = 500;
time = -1:1/srate:1;
frequency = 5;
amplitude = 3;
phase = pi/2;
sinewave = amplitude*sin(2*pi*frequency*time+phase);
figure
subplot(2,1,1)
plot(time,sinewave,'b', LineWidth=3)
hold on

% Generate complex sine wave
complex_sin = amplitude * exp(1i * (2 * pi * frequency * time + phase));
plot(time,complex_sin,'g', LineWidth=3)
xlabel('Time(sec)'), ylabel('Amplitude')
title('Complex sine wave')
legend({'Real';'Imaginary'})


% Generate complex sine wave in 3D
subplot(2,1,2)
plot3(time,real(complex_sin),imag(complex_sin), 'k')
xlabel('Time(Sec)'), ylabel('Real'), zlabel('Imaginary')
axis square
rotate3d on