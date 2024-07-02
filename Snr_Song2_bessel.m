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

Df = Fs / length(data);
f = (-Fs/2:Df:Fs/2-Df);

f_L_range = 3000:1000:10000;
f_H_range = 13000:1000:22000;
SNR_values = zeros(length(f_L_range), length(f_H_range));

for i = 1:length(f_L_range)
    for j = 1:length(f_H_range)
        f_L = f_L_range(i);
        f_H = f_H_range(j);

        %% Bessel bandpass filter
        Order = 4;
        [b_bes, a_bes] = besself(Order/2, 2*pi*[f_L, f_H], 'bandpass');

        %% Bilinear transformation
        [Bz_bes, Az_bes] = bilinear(b_bes, a_bes, Fs);

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

        %% Compute the error signal
        e = x_out_final - x_ref2;

        %% Compute the SNR
        Ps = mean(x_ref2.^2);
        Pn = mean(e.^2); 
        SNR_values(i, j) = 10 * log10(Ps / Pn);  % SNR in dB
    end
end

%% Plot the SNR VS bandpass filter cutoff frequencies 

[X, Y] = meshgrid(f_H_range, f_L_range);
Z = SNR_values;
contour(X, Y, Z,'ShowText','On');
xlabel('f_H (Hz)');
ylabel('f_L (Hz)');
grid on;
title('SNR VS Bandpass Filter Cutoff Frequencies');
