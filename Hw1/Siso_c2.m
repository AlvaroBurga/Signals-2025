%% H1 Alvaro Burga
clc;
clear;
close all;
%% SISO d - function
function [x, x_est,x_est_hReal] = x_identification_c2(L , SNR)

    % First we calculate h in the same way we did it on a.
    N = 200;
    M = L + N - 1;
    variance_w = 1 / SNR;
    theta = 2*pi*rand(M,1) - pi;
    frec = pi/4;
    t= linspace(1,M,M)';
    
    w = sqrt(variance_w*2) * cos(frec*t + theta);

    r = xcorr(w, 'biased'); % To get the autocorrelation of the noise (equal to the covariance as mean is zero)
    rpos = r(M : 2*M-1); % We take from Rw0 to RwM
    Cw = toeplitz(rpos); % We create the noise covariance matrix

    x = randn(N, 1);
    y = zeros(M, 1);
    X = convmtx(x, L);
    
    y(1) = x(1) + w(1);
    for n = 2:N
        y(n) = x(n) + 0.9*y(n-1) + w(n);
    end
    for n = N+1:M
        y(n) = 0.9*y(n-1) + w(n);
    end
    
    h = inv(X' * X)* X' * y;

    % Then we truncate to respect the constrain of len(y) = 200.
    y_truncated = y(1:N);
    H = convmtx(h, N);
    H_truncated = H(1:N, 1:N);
    Cw_truncated = Cw(1:N, 1:N);
    
    % How x is random, we use a bayesian estimator
    x_est = H_truncated' * inv(H_truncated *H_truncated' + Cw_truncated)*y_truncated;

    %Real h
    b = 1;            % MA part
    a = [1 -0.9];     % AR part
    delta = [1; zeros(L-1,1)];  % unit impulse
    h_real = filter(b, a, delta);


    % Then we truncate to respect the constrain of len(y) = 200.
    H = convmtx(h_real, N);
    I = eye(M);
    
    % How x is random, we use a bayesian estimator
    x_est_hReal = H' * inv(H *H' + variance_w*I)*y;
end

%% SISO d results

SNR_dB_values = [-10, -5, 0, 5, 10, 15, 20, 30, 40];  % SNR in dB

SNR_values = 10.^(SNR_dB_values / 10);
L_values = [2,5,10,20,50,100];

% Preallocate MSE matrix
MC = 100;  % number of Monte Carlo runs
MSE = zeros(length(L_values), length(SNR_values));
MSE_hReal = zeros(length(L_values), length(SNR_values));


for i = 1:length(L_values)
    for j = 1:length(SNR_values)
        mse_temp = zeros(1, MC);  % store MSE for each run
        mse_temp_hReal = zeros(1, MC);
        for k = 1:MC
            [x, x_est,x_est_hReal] = x_identification_c2(L_values(i), SNR_values(j));
            mse_temp(k) = mean((x - x_est).^2);
            mse_temp_hReal(k) = mean((x - x_est_hReal).^2);
        end
        MSE(i,j) = mean(mse_temp);  % average over Monte Carlo runs
        MSE_hReal(i,j) = mean(mse_temp_hReal);  % average over Monte Carlo runs
    end
end
MSE_dB = 10 * log10(MSE);
MSE_dB_hReal = 10 * log10(MSE);
%% --- Plot 1: MSE (dB) vs SNR (dB) for different L ---
figure;
hold on;
for i = 1:length(L_values)
    plot(SNR_dB_values, MSE_dB(i,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', ['L = ' num2str(L_values(i))]);
end
xlabel('SNR (dB)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR for different L');
legend('Location', 'northeast');
grid on;
hold off;

% --- Plot 2: MSE (dB) vs L for different SNR (dB) ---
figure;
hold on;
for j = 1:length(SNR_values)
    plot(L_values, MSE_dB(:,j), '-s', 'LineWidth', 1.5, ...
        'DisplayName', ['SNR = ' num2str(SNR_dB_values(j)) ' dB']);
end
xlabel('L (filter length)');
ylabel('Mean Squared Error (dB)');
title('MSE vs L for different SNR');
legend('Location', 'northeast');
grid on;
hold off;

%% --- Plot 2: With h real
figure;
hold on;
for i = 1:length(L_values)
    plot(SNR_dB_values, MSE_dB_hReal(i,:), '-o', 'LineWidth', 1.5, ...
        'DisplayName', ['L = ' num2str(L_values(i))]);
end
xlabel('SNR (dB)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR for different L');
legend('Location', 'northeast');
grid on;
hold off;

% --- Plot 2: MSE (dB) vs L for different SNR (dB) ---
figure;
hold on;
for j = 1:length(SNR_values)
    plot(L_values, MSE_dB_hReal(:,j), '-s', 'LineWidth', 1.5, ...
        'DisplayName', ['SNR = ' num2str(SNR_dB_values(j)) ' dB']);
end
xlabel('L (filter length)');
ylabel('Mean Squared Error (dB)');
title('MSE vs L for different SNR');
legend('Location', 'northeast');
grid on;
hold off;
