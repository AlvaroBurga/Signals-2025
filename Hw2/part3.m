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
[s_full_int, Fs] = audioread('star.wav'); %Reading interference the audio

t_start = 90;   % start time in seconds
t_end   = 94;   % end time in second


i1 = floor(t_start * Fs) + 1;   % starting sample index
i2 = floor(t_end   * Fs);       % ending sample index
s = s_full(i1:i2,1) + s_full(i1:i2,2);
s_int = s_full_int(i1:i2,1) + s_full_int(i1:i2,2);

M = i2 - i1 + 1; %Number of samples evaluated
t = (0:M-1)' / Fs; %Time vector

step = (0:N-1) - (N-1)/2;
x = step*d/(N-1);

A = L; %Long radious of the elipse
B = L/2; % Small radious of the elipse
x_source = A*cos(w0*t); %Position x of the source
y_source = B*sin(w0*t); %Position y of the source
x_int = -2*L * ones(M, 1);
y_int = 2*L * ones(M, 1);

mu = 0.1; %Relative step

points = 100; %Number of points to analyse the ToD
win_size = ceil(M/points); %Window size of the analisis
%% --- Part 2.2.e --- in DToD with Montecarlo
resolution = 20;
SNR_dB_values = [0, 10, 20, 30];
N_mc = 10; % Number of Montecarlo iterations

% Accumulator for Montecarlo
mse_mc = zeros(points, length(SNR_dB_values), N_mc);

for mc = 1:N_mc
    fprintf('Montecarlo iteration: %d/%d\n', mc, N_mc);

    estimated_signal = zeros(M, length(SNR_dB_values));
    filtered_estimated_signal = zeros(M, length(SNR_dB_values));

    signal_recovered = cell(1, N);
    for i = 1:N
        signal_recovered{i} = zeros(M, length(SNR_dB_values));
    end

    for idx = 1:length(SNR_dB_values)
        SNR = SNR_dB_values(idx);
        mics = cell(1, N);
        ToDs = cell(1, N);

        % Generating the rx signals for each mic with the interference
        for mic = 1:N
            [mics{mic}, ToDs{mic}] = rx_generator(SNR, Fs, x(mic), x_source, y_source, vp, s_full, i1, i2);
            if (x(mic) <= 0)
                interference = rx_generator(SNR, Fs, x(mic), x_int, y_int, vp, s_full_int*15, i1, i2);
                mics{mic} = mics{mic} + interference;
            end
        end

        % Calculate the delay for each pair
        mid = ceil(N/2);
        for j = mid:mid
            for k = 1:N
                win = floor(win_size/2);
                for win_i = win:win_size:M-win
                    win = floor(win_size/2);
                    a = max(1, win_i - win);
                    b = min(length(mics{k}), win_i + win);
                    t0 = win_i / Fs;
                    delta_tau_tmp = tau_estimator(t0, mics{j}, mics{k}, Fs, win_size, resolution);
                    signal_recovered{k}(a:b, idx) = recover_signal(mics{j}, delta_tau_tmp, win_i, win_size);
                end
            end
        end

        for i = 1:N
            estimated_signal(:, idx) = estimated_signal(:, idx) + signal_recovered{i}(:, idx);
        end
         estimated_signal(:, idx) = estimated_signal(:, idx) / N;


         %Applying adaptive filter to remove interference
         af_window = win_size*2;
         af = zeros(af_window+1,1);
         win = floor(af_window/2);
         for n = win+1:M-win
             filtered_estimated_signal(n,idx) = af' * estimated_signal(n-win:n+win, idx);
             error = s(n) - filtered_estimated_signal(n,idx);
             trR = estimated_signal(n-win:n+win,idx)'*estimated_signal(n-win:n+win,idx);
             af = af + mu*error*estimated_signal(n-win:n+win,idx)/trR;
         end
         

    end

   
    % Compute MSE blocks for this Montecarlo run
    for idx = 1:length(SNR_dB_values)
        mse_raw = 10*log10((filtered_estimated_signal(:, idx) - s).^2);
        for p = 1:points
            idx_start = (p-1)*win_size + 1;
            idx_end   = min(p*win_size, M);
            mse_mc(p, idx, mc) = mean(mse_raw(idx_start:idx_end));
        end
    end
end

% Average over Montecarlo iterations
mse_avg = mean(mse_mc, 3);


%% --- Part 2.2.c --- MSE
% Time blocks
t_blocks = zeros(points, 1);
for p = 1:points
    idx_start = (p-1)*win_size + 1;
    idx_end   = min(p*win_size, M);
    t_blocks(p) = mean(t(idx_start:idx_end));
end

figure;
hold on;
colors = lines(length(SNR_dB_values));
for idx = 1:length(SNR_dB_values)
    plot(t_blocks, mse_avg(:, idx), '-o', 'LineWidth', 1.5, ...
        'Color', colors(idx,:), ...
        'DisplayName', sprintf('SNR = %d dB', SNR_dB_values(idx)));
end
hold off;
xlabel('Time (s)')
ylabel('Average MSE (dB)')
title('MSE vs Time')
legend show
grid on
%% Sound
r = [estimated_signal(:,4) estimated_signal(:,4)];
sound(r, Fs);