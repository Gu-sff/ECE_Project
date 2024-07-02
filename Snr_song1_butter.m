clear all;
close all;
clc;

%% General parameters
Fs = 44100; % Sampling frequency
load('299097.mat');
x_in = xtot;
x_in = x_in(:,1);

[data, Fs_ref] = audioread('Europe - The Final Countdown.mp3');
x_ref = data(:,1);

% Make sure both signals have the same length
min_length = min(length(x_in), length(x_ref));
x_in = x_in(1:min_length);
x_ref = x_ref(1:min_length);

% Compute SNR in db without filtering
e_o = x_in - x_ref;
Ps_o = mean(x_in.^2);
Pn_o = mean(e_o.^2); 
SNR_o = 10 * log10(Ps_o / Pn_o);

fc_range = 1000:500:21000;
SNR_values = zeros(size(fc_range));

for i = 1:length(fc_range)
    fc = fc_range(i);
    Wo = 2 * pi * fc; 

    %% Butterworth filter (Laplace)
    Order = 4;
    [b_butter, a_butter] = butter(Order,Wo,'s');

    %% Bilinear transformation
    [Bz_butter, Az_butter] = bilinear(b_butter, a_butter, Fs);

    %% Apply the digital filter to the signal
    x_out_butter = filter(Bz_butter, Az_butter, x_in);

    %% Realignment
    delay = finddelay(x_ref, x_out_butter);
    x_out_butter = circshift(x_out_butter, -delay);

    %% Compute the error signal
    e = x_out_butter - x_ref;

    %% Compute the SNR
    Ps = mean(x_ref.^2);  % Signal power
    Pn = mean(e.^2);  % Noise power
    SNR_values(i) = 10 * log10(Ps / Pn);  % SNR in dB
end

%% Plot SNR versus filter bandwidth
figure;
plot(fc_range, SNR_values, '-o', 'DisplayName', 'Filtered');
hold on;
yline(SNR_o, 'r-', 'DisplayName', 'Unfiltered');
grid on;
xlabel('Cutoff Frequency (Hz)');
ylabel('SNR (dB)');
title('SNR vs. Filter Bandwidth');
legend;
hold off;