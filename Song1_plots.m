clear all;
close all;
clc;

%% General parameters
Fs = 44100; % Sampling frequency
load('299097.mat');
data = xtot;
[data_ref, Fs_ref] = audioread('Europe - The Final Countdown.mp3');
x_ref = data_ref(:,1);
Ts=1/Fs;
t=(0:Ts:(length(data)-1)*Ts);

% Make sure both signals have the same length
x_in = data(:,1);
min_length = min(length(x_in), length(x_ref));
x_in = x_in(1:min_length);
x_ref = x_ref(1:min_length);

%% Plot of original signal spectrum
Df=Fs/length(data);
f=(-Fs/2:Df:Fs/2-Df);
F=fftshift(fft(x_in));
figure;
plot(f,10*log10(abs(F).^2));
hold on;
grid on;
xlabel('Hz')
ylabel('[dB[A.U.]')
title('Plot of original signal spectrum')
grid on
%% Bessel filter (Laplace)
Order=4;
fc_bes=7500; % optimal cut off frequency 
Wo=2*pi*fc_bes;
[b_bes,a_bes] = besself(Order,Wo);

%% Analog Bessel filter Fourier transform
s = 1j*2*pi*f;
Ha_bes = polyval(b_bes,s)./polyval(a_bes,s);
Ha_bes = Ha_bes/max(Ha_bes); 

%% Bilinear transformation
[Bz_bes,Az_bes]=bilinear(b_bes,a_bes,Fs);

%% Digital Bessel filter discrete fourier transform
z=exp(1j*2*pi*f/Fs);
Hz_bes = polyval(Bz_bes,z)./polyval(Az_bes,z);
Hz_bes = Hz_bes/max(Hz_bes);

%% Apply the digital filter to the signal
x_out_bes = filter(Bz_bes, Az_bes, x_in);

%% Plot the Analog Bessel filter Fourier transform
figure;
plot(f, 20 * log10(abs(Ha_bes)));
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold on 
%% Plot the Digital Bessel filter discrete Fourier transform
plot(f, 20 * log10(abs(Hz_bes)));
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('Analog Bessel Filter Fourier Transform (Ha_bes)','Digital Bessel Filter Discrete Fourier Transform (Hz_bes)')
title('Transfer functions for Bessel filter');
grid on;

%% Plot the filtered signal in the frequency domain
F_out_bes = fftshift(fft(x_out_bes));
figure;
plot(f, 10 * log10(abs(F_out_bes).^2));
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('Frequency Domain: Filtered Signal(Bessel)');

%% Butterworth filter (Laplace)
Order=4; 
fc_butter=6500; % optimal cut off frequency 
Wn = 2*pi*fc_butter;
[b_butter,a_butter] = butter(Order,Wn,'s');

%% Analog Butterworth filter Fourier transform
s = 1j*2*pi*f;
Ha_butter = polyval(b_butter,s)./polyval(a_butter,s);
Ha_butter = Ha_butter/max(Ha_butter); 

%% Bilinear transformation
[Bz_butter,Az_butter]=bilinear(b_butter,a_butter,Fs);

%% Digital Butterworth filter discrete fourier transform
z=exp(1j*2*pi*f/Fs);
Hz_butter = polyval(Bz_butter,z)./polyval(Az_butter,z);
Hz_butter = Hz_butter/max(Hz_butter);

%% Apply the digital filter to the signal
x_out_butter = filter(Bz_butter, Az_butter, x_in);

%% Plot the Analog Butterworth filter Fourier transform
figure;
plot(f, 20 * log10(abs(Ha_butter)));
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold on

%% Plot the Digital Butterworth filter discrete Fourier transform
plot(f, 20 * log10(abs(Hz_butter)));
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('Frequency Response of Analog Butterworth Filter(Ha_butter)','Frequency Response of Digital Butterworth Filter(Hz_butter)')
title('Transfer functions for Butterworth filter');

%% Frequency Domain Plot of Filtered Signal using FFT
F_out_butter = fftshift(fft(x_out_butter));
figure;
plot(f, 10 * log10(abs(F_out_butter).^2));
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('Frequency Domain: Filtered Signal(Butterworth)');

%%
figure;
plot(t,x_ref);
hold on;
plot(t,x_out_bes);
grid on;
xlabel('time [s]');
ylabel('amplitude');
legend('Reference Signal', 'Error Signal with filter')
title('In time domain: Reference Signal Vs. Error Signal(Bessel low-pass filter applied)')

%%
figure;
plot(t,x_ref);
hold on;
plot(t,x_out_butter);
grid on;
xlabel('time [s]');
ylabel('amplitude');
legend('Reference Signal', 'Error Signal with filter')
title('In time domain: Reference Signal Vs. Error Signal(Butterworth low-pass filter applied')

