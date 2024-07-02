clear all;
close all;
clc;

%% General parameters
Fs = 44100; 
load('299097.mat');
data = xtot(:, 1);
[data_song1, fs] = audioread('Europe - The Final Countdown.mp3');
x_ref1 = data_song1(:, 1);
[data_song2, fs2] = audioread('A-Ha_-_Take_On_Me.mp3');
x_ref2 = data_song2(:, 1);

% Make sure both signals have the same length
min_length = min(length(data), length(x_ref1));
data = data(1:min_length);
x_ref1 = x_ref1(1:min_length);

% Make sure both signals have the same length
min_length = min(length(data), length(x_ref2));
data = data(1:min_length);
x_ref2 = x_ref2(1:min_length);

Df=Fs/length(data);
f=(-Fs/2:Df:Fs/2-Df);
F=fftshift(fft(data));
F1=fftshift(fft(x_ref1));
F2=fftshift(fft(x_ref2));

figure;
plot(f,10*log10(abs(F).^2));
hold on;
grid on;
xlabel('Hz')
ylabel('[dB[A.U.]')
title('Plot of error signal spectrum')
grid on

figure;
plot(f,10*log10(abs(F1).^2));
hold on;
grid on;
xlabel('Hz')
ylabel('[dB[A.U.]')
title('Plot of reference signal(song1) spectrum')
grid on

figure;
plot(f,10*log10(abs(F2).^2));
hold on;
grid on;
xlabel('Hz')
ylabel('[dB[A.U.]')
title('Plot of reference signal(song2) spectrum')
grid on



