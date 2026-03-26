%% H2 Alvaro Burga
clc;
clear;
close all;

%% Setting
N = 15; % Number of antennas
vp = 343; % Propagation speed
d = 0.18; % Distance of the edge antennas
L=2; %Long radious of the elipsis
w0 = pi/2; %Angular frecuency of the rotation

[s_full, Fs] = audioread('lavanda.wav'); %Reading the audio
t_start = 50;   % start time in seconds
t_end   = 54;   % end time in seconds

i1 = floor(t_start * Fs) + 1;   % starting sample index
i2 = floor(t_end   * Fs);       % end,ing sample index
s = s_full(i1:i2,1) + s_full(i1:i2,2);

M = i2 - i1 + 1; %Number of samples evaluated
t = (0:M-1)' / Fs; %Time vector

step = (0:N-1) - (N-1)/2;
x = step*d/(N-1);

A = L; %Long radious of the elipse
B = L/2; % Small radious of the elipse
x_source = A*cos(w0*t); %Position x of the source
y_source = B*sin(w0*t); %Position y of the source 

points = 100; %Number of points to analyse the ToD
win_size = ceil(M/points); %Window size of the analisis
%% --- Part 2.2.a --- in DToD
resolution = 20;
SNR_dB_values = [-10, 0, 10, 20, 30];
R = zeros(M,length(SNR_dB_values), N);
MC = 20; %Montecarlo simulation
MSE_DToD= zeros(M,length(SNR_dB_values));
delta_tau = zeros(M,length(SNR_dB_values));
signal_recovered = zeros(M,length(SNR_dB_values), N);

for j = 1:N
    for k = (j+1):N
        for idx = 1:length(SNR_dB_values)
        
            SNR = SNR_dB_values(idx);
        
            delta_tau_tmp= zeros(M,MC);
            MSE_DToD_tmp= zeros(M,MC);
            signal_recovered_tmp_j = zeros(M,MC);
            signal_recovered_tmp_k = zeros(M,MC);
            for m = 1:MC
                [mic_left,ToD_left] = rx_generator(SNR,Fs,x(j),x_source,y_source,vp,s_full,i1,i2);
                [mic_right,ToD_right] = rx_generator(SNR,Fs,x(k),x_source,y_source,vp,s_full,i1,i2);
        
                K1 = floor(win_size/2) + 1 ;
                K2 = K1;
        
                for i = K1:win_size:M-K2
        
                    win = floor(win_size/2);
                    a = max(1, i - win);
                    b = min(length(mic_left), i + win);
        
                    t0 = i /Fs;
                    delta_tau_tmp(i,m) = delta_tau_tmp(i,m) + tau_estimator(t0, mic_left, mic_right, Fs, win_size, resolution);
                    real_tau_tmp = ToD_right(i) -  ToD_left(i);
                    MSE_DToD_tmp(i,m) = MSE_DToD_tmp(i,m) + (real_tau_tmp - delta_tau_tmp(i,m))^2;

                    mp = 1;
                    %Use this to compensate the attenuation
                    %mp = meanPower(mic_left, mic_right, win_size, i);

                    signal_recovered_tmp_j(i-win_size/2 : i + win_size/2,m) = recover_signal(mic_left, delta_tau_tmp(i,m), i, win_size, mp);
                    signal_recovered_tmp_k(i-win_size/2 : i + win_size/2,m) = recover_signal(mic_right, -delta_tau_tmp(i,m), i, win_size, mp);

                end
            end
            signal_recovered(:, idx,j) = signal_recovered(:, idx,j) + mean(signal_recovered_tmp_j, 2);
            signal_recovered(:, idx,k) = signal_recovered(:, idx,k) + mean(signal_recovered_tmp_k, 2);
            delta_tau(:,idx) = mean(delta_tau_tmp, 2);
            MSE_DToD(:, idx) = mean(MSE_DToD_tmp, 2);
            real_tau = ToD_right - ToD_left;
        end
    end
end

signal_recovered = signal_recovered/(N-1);
%% --- Plot 2.1.b: time (s) vs time delay (s) for different SNR(dB) ---
MSE_dB = 10*log(MSE_DToD);
numSNR = length(SNR_dB_values);

for i = 1:numSNR

    subplot(numSNR, 1, i)  
    plot(t, real_tau, 'k', 'LineWidth', 1.5); hold on
    plot(t, delta_tau(:,i), 'LineWidth', 1.5);

    title(['SNR = ' num2str(SNR_dB_values(i)) ' dB'])
    ylabel('\tau (s)')
    grid on

    if i == numSNR
        xlabel('t (s)')
    end

    legend('Real delay','Estimated delay','Location','best')
end

%% 2.1.b. To hear the recoverd signal of the max power
R= sum(signal_recovered,3)/N;
mse = 10*log((R - s).^2);
mean_mse = mean(R, 1);
figure;

plot(SNR_dB_values, mean_mse, '-o', 'LineWidth', 1.5)
xlabel('SNR (dB)')
ylabel('Average MSE')
title('average')
grid on

%% Sound
r = [R(:,5) R(:,5)];
sound(r, Fs);