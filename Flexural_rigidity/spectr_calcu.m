function [ psdx,freq,relaxationtime] = spectr_calcu( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%  x=x*1e-9
Fs = 19.7;
t = 0:1/Fs:1-1/Fs;


N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;
     [psdx,freq] = pwelch(x,[],[],[],2);
figure (467); loglog(freq(3:end),(psdx(3:end)))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency ')
nlinearfitting( freq(1:end),(psdx(1:end)) )
% [fitresult, gof] = fourier_spectrum_fit(freq(3:end),(psdx(3:end)))
  [tr, tr_std1, A, pdf] = fit_spectrum_nonli(freq(2:12),(psdx(2:12)))
relaxationtime=tr;
end

