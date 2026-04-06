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
%% --- Part 2.2.c --- in DToD
resolution = 20;
SNR_dB_values = [-10, 0, 10, 20, 30, 40];
estimated_signal= zeros(M,length(SNR_dB_values));

signal_recovered = cell(1, N); %The Signal recovered for each mic
for i = 1:N
    signal_recovered{i} = zeros(M, length(SNR_dB_values)); %The rx signal for each SNR
end


for idx = 1:length(SNR_dB_values) %For each SNR value defined

    SNR = SNR_dB_values(idx);
    mics = cell(1, N); %The recieved signals
    ToDs = cell(1, N); %The delay from the source for each mic

    %Generating the rx signals for each mic
    for mic = 1:N
        [mics{mic}, ToDs{mic}] = rx_generator(SNR,Fs,x(mic),x_source,y_source,vp,s_full,i1,i2);
    end

    %Calculate the dalay for each pair
    mid = ceil(N/2);
    for j = mid:mid
        for k = 1:N

            win = floor(win_size/2);
            for win_i =win:win_size:M-win %For each time instant

                win = floor(win_size/2);
                a = max(1, win_i - win); %First point of the grouping
                b = min(length(mics{k}), win_i + win); %last point of the grouping

                t0 = win_i /Fs; %The time asociated with the group
                delta_tau_tmp = tau_estimator(t0, mics{j}, mics{k}, Fs, win_size, resolution);
                signal_recovered{k}(a : b,idx) = recover_signal(mics{j}, delta_tau_tmp, win_i, win_size);

            end

        end
    end

    for i = 1:N
        estimated_signal(:,idx) = estimated_signal(:,idx) + signal_recovered{i}(:,idx);
    end
end

estimated_signal = estimated_signal/N;



%% --- Part 2.2.c --- MSE
mse = 10*log((estimated_signal - s).^2);
mean_mse = mean(mse, 1);
figure;

plot(SNR_dB_values, mean_mse, '-o', 'LineWidth', 1.5)
xlabel('SNR (dB)')
ylabel('Average MSE')
title('MSE vs SNR')
grid on

%%  --- Part 2.2.d --- in DToD
resolution = 20;
SNR_dB_values = [-10, 0, 10, 20, 30, 40];
estimated_signal= zeros(M,length(SNR_dB_values));

signal_recovered = cell(1, N); %The Signal recovered for each mic
for i = 1:N
    signal_recovered{i} = zeros(M, length(SNR_dB_values)); %The rx signal for each SNR
end


for idx = 1:length(SNR_dB_values) %For each SNR value defined

    SNR = SNR_dB_values(idx);
    mics = cell(1, N); %The recieved signals
    ToDs = cell(1, N); %The delay from the source for each mic

    %Generating the rx signals for each mic

    %Noise
    SNR_linear = 10^(SNR/10); 
    var_w = 1/SNR_linear;
    Rw = zeros(N,N);
    for i = 1:N
        for j =1:N
            if(i == j)
                Rw(i,j) = var_w;
            elseif(abs(i-j) == 1)
                Rw(i,j)= var_w*0.2;
            else
                Rw(i,j) = 0;
            end
        end
    end

    [eigenvectors, eigenvalues] = eig(Rw);

    w = eigenvectors*sqrt(eigenvalues)*randn(N,M); %Correlated noise

    for mic = 1:N
        [mics{mic}, ToDs{mic}] = rx_generator_no_noise(Fs,x(mic),x_source,y_source,vp,s_full,i1,i2);
        mics{mic} = mics{mic} + w(mic,:)';
    end

    %Calculate the dalay for each pair
    mid = ceil(N/2);
    for j = 1:1
        for k = 1:N

            win = floor(win_size/2);
            for win_i =win:win_size:M-win %For each time instant

                win = floor(win_size/2);
                a = max(1, win_i - win); %First point of the grouping
                b = min(length(mics{k}), win_i + win); %last point of the grouping

                t0 = win_i /Fs; %The time asociated with the group
                delta_tau_tmp = tau_estimator(t0, mics{j}, mics{k}, Fs, win_size, resolution);
                signal_recovered{k}(a : b,idx) = recover_signal(mics{j}, delta_tau_tmp, win_i, win_size);

            end

        end
    end

    gamma = gamma_weights(N, 2);
    for i = 1:N
        estimated_signal(:,idx) = estimated_signal(:,idx) + signal_recovered{i}(:,idx)*gamma(i);
    end
end


%% --- Part 2.2.d --- MSE
mse = 10*log((estimated_signal - s).^2);
mean_mse = mean(mse, 1);
figure;

plot(SNR_dB_values, mean_mse, '-o', 'LineWidth', 1.5)
xlabel('SNR (dB)')
ylabel('Average MSE')
title('MSE vs SNR')
grid on

%% Sound
r = [estimated_signal(:,5) estimated_signal(:,5)];
sound(r, Fs);