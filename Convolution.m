load v1_laminar.mat
signal = mean(csd(7,:,:),3);
subplot(4,1,1)
plot(timevec,signal,'k')
set(gca,"XLim",[-.5,1.5])

% Create Gaussian
time = -1:1/srate:1;
fwhm = .5/100;                                         % width of the gausssian in seconds
gaussian_kernel = exp(-time.^2./(fwhm^2));
gaussian_kernel = gaussian_kernel./sum(gaussian_kernel);
subplot(4,1,2)
plot(time,gaussian_kernel)
set(gca,"XLim",[-.5,1.5])

% Frequency domain, both kernel and signal
convlength = length(gaussian_kernel) + length(signal) -1;
signal_specteral = fft(signal,convlength);
kernel_specteral = fft(gaussian_kernel,convlength);
trim = floor(length(gaussian_kernel)/2); 
convresx = ifft(signal_specteral .* kernel_specteral); % Elementwise Multiplication, then get back to time-domain
convresx = convresx(trim+1:end-trim+1);                % The convolved signal its longer than the origninal signal
subplot(4,1,3)
plot(timevec,convresx,'k')
set(gca,"XLim",[-.5,1.5])
subplot(4,1,4)
plot(timevec,signal,'b')
set(gca,"XLim",[-.5,1.5])
hold on
plot(timevec,convresx,'r')
set(gca,"XLim",[-.5,1.5])
legend('Original ERP', 'Gaussian convolved')