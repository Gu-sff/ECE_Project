clear all;
close all;
clc;

%% General parameters
Fs = 44100; 
ID = 299097;
fc_ca = 9500 + floor(ID/1000) + (ID - floor(ID/1000)*1000);
load('299097.mat');
data = xtot(:, 1);
[data_song2, fs] = audioread('A-Ha_-_Take_On_Me.mp3');
x_ref2 = data_song2(:, 1);

% Make sure both signals have the same length
min_length = min(length(data), length(x_ref2));
data = data(1:min_length);
x_ref2 = x_ref2(1:min_length);

Ts = 1/Fs;
t = (0:min_length-1) * Ts;
Df = Fs/length(data);
f = (-Fs/2:Df:Fs/2-Df);
%% Bessel bandpass filter
Order = 4;
f_L = 7000; % Optimal f_L
f_H = 22000; % Optimal f_H
[b_bes, a_bes] = besself(Order/2, 2*pi*[f_L, f_H], 'bandpass');

%% Analog Bessel filter Fourier transform
s = 1j*2*pi*f;
Ha_bes = polyval(b_bes,s)./polyval(a_bes,s);
Ha_bes = Ha_bes/max(Ha_bes); 

%% Bilinear transformation
[Bz_bes, Az_bes] = bilinear(b_bes, a_bes, Fs);

%% Digital Bessel filter discrete fourier transform
z=exp(1j*2*pi*f/Fs);
Hz_bes = polyval(Bz_bes,z)./polyval(Az_bes,z);
Hz_bes = Hz_bes/max(Hz_bes);

%% Apply the digital filter to the signal
x_out_bes = filter(Bz_bes, Az_bes, data);

%% Demodulate the signal
x_out_bes = x_out_bes*2.*cos(2 * pi * fc_ca * t');

%% Low-pass filter after demodulation
fc_bes = 7500;  % Best cut off frequency for the low-pass filter
Wo = 2 * pi * fc_bes;
[b, a] = besself(Order, Wo);
[Bz_lp, Az_lp] = bilinear(b, a, Fs);
x_out_final = filter(Bz_lp, Az_lp, x_out_bes);

%% Realignment
delay = finddelay(x_ref2, x_out_final);
x_out_final = circshift(x_out_final, -delay);

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
title('Transfer functions for optimal Bessel filter');
grid on;

%% Plot reference signal and error signal in optimal bandwidth(Bessel) in time domain
figure;
plot(t,x_ref2);
grid on;
xlabel('time [s]');
ylabel('amplitude [A.U.]');
hold on;
plot(t,x_out_final);
grid on;
xlabel('time [s]');
ylabel('amplitude [A.U.]');
legend('Reference signal in time domain','Error signal in optimal bandwidth in time domain')
title('Reference signal vs the error signal in optimal conditions(Bessel filter applied)')


%% Butterworth bandpass filter
Order = 4;
f_L_butter = 8000;
f_H_butter = 22000;
[b_butter, a_butter] = butter(Order/2, 2*pi*[f_L_butter, f_H_butter], 'bandpass','s');

%% Analog Butterworth filter Fourier transform
s = 1j*2*pi*f;
Ha_butter = polyval(b_butter,s)./polyval(a_butter,s);
Ha_butter = Ha_butter/max(Ha_butter); 

%% Bilinear transformation
[Bz_butter, Az_butter] = bilinear(b_butter, a_butter, Fs);

%% Digital Butterworth filter discrete fourier transform
z=exp(1j*2*pi*f/Fs);
Hz_butter = polyval(Bz_butter,z)./polyval(Az_butter,z);
Hz_butter = Hz_butter/max(Hz_butter);

%% Apply the digital filter to the signal
x_out_butter = filter(Bz_butter, Az_butter, data);

%% Demodulate the signal
x_out_butter = x_out_butter*2.*cos(2 * pi * fc_ca * t');

%% Low-pass filter after demodulation
fc_butter = 6500; % optimal cut off frequency 
Wn = 2*pi*fc_butter;
[b_butter,a_butter] = butter(Order,Wn,'s');
[Bz_butter, Az_butter] = bilinear(b_butter, a_butter, Fs);
x_out_final_butter = filter(Bz_butter, Az_butter, x_out_butter);

%% Realignment
delay = finddelay(x_ref2, x_out_final_butter);
x_out_final_butter = circshift(x_out_final_butter, -delay);

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
title('Transfer functions for optimal Butterworth filter');

%% Plot reference signal and error signal in optimal bandwidth(Butterworth) in time domain
figure;
plot(t,x_ref2);
grid on;
xlabel('time [s]');
ylabel('amplitude [A.U.]');
hold on;
plot(t,x_out_final_butter);
grid on;
xlabel('time [s]');
ylabel('amplitude [A.U.]');
legend('Reference signal in time domain','Error signal in optimal bandwidth in time domain')
title('Reference signal vs the error signal in optimal conditions(Butterworth filter applied)')

