%% H2 Alvaro Burga
clc;
clear;
close all;

%% Setting
N = 2; % Number of antennas
vp = 343; % Propagation speed
d = 0.18; % Distance of the antennas
L=2; %Long radious of the elipsis
w0 = pi/2; %Angular frecuency of the rotation

[s_full, Fs] = audioread('lavanda.wav'); %Reading the audio
t_start = 50;   % start time in seconds
t_end   = 54;   % end time in seconds

i1 = floor(t_start * Fs) + 1;   % starting sample index
i2 = floor(t_end   * Fs);       % ending sample index
s = s_full(i1:i2,1) + s_full(i1:i2,2);

M = i2 - i1 + 1; %Number of samples evaluated
t = (0:M-1)' / Fs; %Time vector

x_left = -d/2; %Position of the left mic
x_right = d/2; %Position of the right mic
A = L; %Long radious of the elipse
B = L/2; % Small radious of the elipse
x_source = A*cos(w0*t); %Position x of the source
y_source = B*sin(w0*t); %Position y of the source 

points = 100; %Number of points to analyse the ToD
win_size = ceil(M/points); %Window size of the analisis

%% Part 2.1.a
resolution = 20; % Number of samples to interpolate
t0_values = [1, 2.5, 2];
SNR_dB_values = [-30, -10, 0, 10, 30, 50, 70];
MC = 300;

tau_est = zeros(length(t0_values), length(SNR_dB_values));
tau_real = zeros(length(t0_values), length(SNR_dB_values));
mse = zeros(length(t0_values), length(SNR_dB_values));
for i = 1:length(t0_values)
    index = floor(t0_values(i)* Fs) +1;
    for j = 1:length(SNR_dB_values)
        tau_temp = zeros(MC,1);
        mse_temp = zeros(MC,1);
        for m = 1:MC
            [mic_left,ToD_left] = rx_generator(SNR_dB_values(j),Fs,x_left,x_source,y_source,vp,s_full, i1, i2);
            [mic_right,ToD_right] = rx_generator(SNR_dB_values(j),Fs,x_right,x_source,y_source, vp,s_full, i1, i2);  
            tau_real(i,j) =  ToD_right(index) - ToD_left(index);
            tau_temp(m) = tau_estimator(t0_values(i), mic_left, mic_right, Fs,win_size, resolution);
            error_time = tau_temp(m) - tau_real(i,j);
            mse_temp(m) = error_time^2;
        end
        mse(i,j) = mean(mse_temp);
        tau_est(i,j) = mean(tau_temp);
    end
end
rmse = sqrt(mse);

%% --- Plot 2.1.a: RMSE (s) vs SNR (dB) for different t(s) ---
figure;
hold on;
for i = 1:length(t0_values)
    plot(SNR_dB_values, rmse(i,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', ['t = ' num2str(t0_values(i))]);
end
xlabel('SNR (dB)');
ylabel('Root Mean Squared Error  (s)');
title('MSE vs SNR for different t0');
legend('Location', 'northeast');
grid on;
hold off;
%% --- Part 2.1.b --- in DToD
resolution = 20;
SNR_dB_values = [-10, 0, 10, 20, 30];
R_right = zeros(M,length(SNR_dB_values));
R_left = zeros(M,length(SNR_dB_values));
MC = 30; %Montecarlo simulation
MSE_DToD= zeros(M,length(SNR_dB_values));
delta_tau = zeros(M,length(SNR_dB_values));
signal_recovered_L = zeros(M,length(SNR_dB_values));
signal_recovered_R = zeros(M,length(SNR_dB_values));


for idx = 1:length(SNR_dB_values)

    SNR = SNR_dB_values(idx);

    delta_tau_tmp= zeros(M,MC);
    MSE_DToD_tmp= zeros(M,MC);
    signal_recovered_tmp_L = zeros(M,MC);
    signal_recovered_tmp_R = zeros(M,MC);

    for m = 1:MC
        [mic_left,ToD_left] = rx_generator(SNR,Fs,x_left,x_source,y_source,vp,s_full,i1,i2);
        [mic_right,ToD_right] = rx_generator(SNR,Fs,x_right,x_source,y_source,vp,s_full,i1,i2);

        win = floor(win_size/2);

        for i = win:win_size:M-win %Working in a window, to not consider all points
            
            a = max(1, i - win);
            b = min(length(mic_left), i + win);

            t0 = i /Fs;
            delta_tau_tmp(i,m) = delta_tau_tmp(i,m) + tau_estimator(t0, mic_left, mic_right, Fs, win_size, resolution);
            real_tau_tmp = ToD_right(i) -  ToD_left(i);
            MSE_DToD_tmp(i,m) = MSE_DToD_tmp(i,m) + (real_tau_tmp - delta_tau_tmp(i,m))^2;

            signal_recovered_tmp_L(a:b,m) = recover_signal(mic_left, delta_tau_tmp(i,m), i, win_size);
            signal_recovered_tmp_R(a:b,m) = recover_signal(mic_right, delta_tau_tmp(i,m)*-1, i, win_size);

        end
    end
    signal_recovered_L(:, idx) = mean(signal_recovered_tmp_L, 2);
    signal_recovered_R(:, idx) = mean(signal_recovered_tmp_R, 2);
    delta_tau(:,idx) = mean(delta_tau_tmp, 2);
    MSE_DToD(:, idx) = mean(MSE_DToD_tmp, 2);
    real_tau = ToD_right - ToD_left;
end

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
mse_L = 10*log((signal_recovered_L - s).^2);
mse_R = 10*log((signal_recovered_R - s).^2);
mean_mse_L = mean(mse_L, 1);
mean_mse_R = mean(mse_R, 1);

figure;

subplot(1,2,1)
plot(SNR_dB_values, mean_mse_L, '-o', 'LineWidth', 1.5)
xlabel('SNR (dB)')
ylabel('Average MSE')
title('Left channel')
grid on

subplot(1,2,2)
plot(SNR_dB_values, mean_mse_R, '-o', 'LineWidth', 1.5)
xlabel('SNR (dB)')
ylabel('Average MSE')
title('Right channel')
grid on


%% Sound
r = [signal_recovered_L(:,1) signal_recovered_R(:,1)];
sound(r, Fs);