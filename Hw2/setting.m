%% H2 Alvaro Burga
clc;
clear;
close all;

%% Setting
N = 2;
vp = 343;
d = 0.18;
L=2;
w0 = pi/2;
maxLagTime = 0.001;

[s_full, Fs] = audioread('lavanda.wav');
t_start = 50;   % start time in seconds
t_end   = 54;   % end time in seconds

i1 = floor(t_start * Fs) + 1;   % starting sample index
i2 = floor(t_end   * Fs);       % ending sample index

M = i2 - i1 + 1;
t = (0:M-1)' / Fs;

x_left = -d/2;
x_right = d/2;
A = L;
B = L/2;
x_source = A*cos(w0*t);
y_source = B*sin(w0*t);

%% Analysis
mic_left= rx_generator(-10,Fs,x_left,x_source,y_source,vp,s_full, i1, i2);

mic_right= rx_generator(-10,Fs,x_right,x_source,y_source,vp,s_full, i1, i2);
figure(1)
plot(t, mic_left);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time Evolution of Left Channel');

Y = fft(mic_left);        % Use left channel for FFT
Y = Y(1:M/2);           % keep positive frequencies only
f = (0:M/2-1) * (Fs/M); % frequency axis

figure(2)
plot(f, abs(Y));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum (15 to 23 s)');
grid on;

figure(3)
plot(f, 20*log10(abs(Y)));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectrum in dB (15 to 23 s)');
grid on;

window = 1024;
noverlap = 512;
nfft = 1024;

figure(4)
spectrogram(mic_left, window, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram (15 to 23 s)');
colorbar;

mic = [mic_left, mic_right];
sound(mic, Fs);

%% Part 2.1
resolution = 0; %Change accordingly
t0_values = [1, 2.5, 2];
SNR_dB_values = [-30, -10, 0, 10, 30, 50, 70];
MC = 20;

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
            tau_temp(m) = tau_estimator(t0_values(i), mic_left, mic_right, Fs,maxLagTime, resolution); 
            mse_temp(m) = (tau_temp(m) - tau_real(i,j))^2;
        end
        mse(i,j) = mean(mse_temp);
        tau_est(i,j) = mean(tau_temp);
    end
end
MSE_dB = 10*log(mse);

%% --- Plot 2.1: MSE (dB) vs SNR (dB) for different t(s) ---
figure;
hold on;
for i = 1:length(t0_values)
    plot(SNR_dB_values, MSE_dB(i,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', ['t = ' num2str(t0_values(i))]);
end
xlabel('SNR (dB)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR for different t0');
legend('Location', 'northeast');
grid on;
hold off;
%% --- Part 2.2 --- in DToD
resolution = 0;
SNR_dB_values = [-30, -10, 0, 10, 30, 50, 70];
source = s_full(i1:i2,1) + s_full(i1:i2,2);
R_right = zeros(M,length(SNR_dB_values));
R_left = zeros(M,length(SNR_dB_values));
MC = 20; %Montecarlo simulation
MSE_DToD= zeros(M,length(SNR_dB_values));
delta_tau = zeros(M,length(SNR_dB_values));
real_tau = zeros(M,length(SNR_dB_values));


for idx = 1:length(SNR_dB_values)

    SNR = SNR_dB_values(idx);

    for m = 1:MC
        [mic_left,ToD_left] = rx_generator(SNR,Fs,x_left,x_source,y_source,vp,s_full,i1,i2);
        [mic_right,ToD_right] = rx_generator(SNR,Fs,x_right,x_source,y_source,vp,s_full,i1,i2);

        K1 = floor((maxLagTime*4) * Fs) + 1;
        K2 = K1;
        delta_tau_tmp= zeros(M,MC);
        MSE_DToD_tmp= zeros(M,MC);
        for i = K1*2:M-2*K2
            t0 = i /Fs;
            delta_tau_tmp(i,m) = tau_estimator(t0, mic_left, mic_right, Fs, maxLagTime, resolution);
            real_tau_tmp = ToD_right(i) -  ToD_left(i);
            MSE_DToD_tmp(i,m) = (real_tau_tmp - delta_tau_tmp(i,m))^2;
        end
    end
    delta_tau(:,idx) = mean(delta_tau_tmp, 2);
    MSE_DToD(:, idx) = mean(MSE_DToD_tmp, 2);
    real_tau = ToD_right - ToD_left;
end

%% --- Plot 2.2: MSE (dB) vs time (s) for different SNR(dB) ---
MSE_dB = 10*log(MSE_DToD);
figure;
hold on;
for i = 1:1%length(SNR_dB_values)
    plot(t, real_tau(:,i), 'LineWidth', 1.5, ...
        'DisplayName', ['Real_SNR = ' num2str(SNR_dB_values(i))]);
        plot(t, delta_tau(:,i), 'LineWidth', 1.5, ...
        'DisplayName', ['Real_SNR = ' num2str(SNR_dB_values(i))]);
end
xlabel('t (s)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR for different SNR');
legend('Location', 'northeast');
grid on;
hold off;


%% --- Part 2.2 adapting the filter in amplitude
SNR_dB_values = [-70 -50 -30 -10 0 10 30 50 70];
mu = 0.4;
MC = 20;

% Resultados
All_R_left  = zeros(M, length(SNR_dB_values));
All_R_right = zeros(M, length(SNR_dB_values));
MSE_left  = zeros(1, length(SNR_dB_values));
MSE_right = zeros(1, length(SNR_dB_values));

source = s_full(i1:i2,1) + s_full(i1:i2,2);

for idx = 1:length(SNR_dB_values)

    SNR = SNR_dB_values(idx);
    R_right = zeros(M,1);
    R_left = zeros(M,1);

    for m = 1:MC
        [mic_left,ToD_left] = rx_generator(SNR,Fs,x_left,x_source,y_source,vp,s_full,i1,i2);
        [mic_right,ToD_right] = rx_generator(SNR,Fs,x_right,x_source,y_source,vp,s_full,i1,i2);

        [R_right_i, R_left_i] = audio_estimator(maxLagTime, M, Fs, source, mu, mic_left, mic_right);

        R_right = R_right + R_right_i;
        R_left = R_left + R_left_i;
    end

    R_right = R_right ./ MC;
    R_left  = R_left  ./ MC;

    All_R_left(:,idx)  = R_left;
    All_R_right(:,idx) = R_right;

    MSE_left(idx)  = mean( (source - R_left).^2 );
    MSE_right(idx) = mean( (source - R_right).^2 );
end

MSE_dB_left = 10*log(MSE_left);
MSE_dB_right = 10*log(MSE_right);
%%

figure(2);
hold on;
plot(SNR_dB_values, MSE_dB_left, '-o', 'LineWidth', 1.5, ...
    'DisplayName', ['MSE Left']);
plot(SNR_dB_values, MSE_dB_right, '-o', 'LineWidth', 1.5, ...
    'DisplayName', ['MSE right']);
xlabel('SNR (dB)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR');
legend('Location', 'northeast');
    
r = [R_right R_left];
sound(r, Fs);


%%
bin_width = 5e-5;

edges = min(tod):bin_width:max(tod);

histogram(tod, edges)