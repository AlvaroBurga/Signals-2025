%% H1 Alvaro Burga
clc;
clear;
close all;

%% Function
function [w, CW] = noise_generator(rho, N)
    L=200;
    variance_w = 1;

    CW = zeros(N,N);
    for j = 1:N
        for i = 1:N
            CW(i,j) = variance_w*rho^(abs(i-j));
        end
    end
    [V,D] = eig(CW);
    w = zeros(L,N);
    for k = 1:L
        z = randn(N,1);
        w(k,:) = V* sqrt(D) * z;
    end
end

%% Result
N_values = linspace(20,200,10);
rho_values = [0, 0.5, 0.9];

MC = 200;  % number of Monte Carlo runs
MSE = zeros(length(rho_values), length(N_values));

for i = 1:length(rho_values)
    for j = 1:length(N_values)
        mse_temp = zeros(1, MC);  % store MSE for each run
        for k = 1:MC
            [w, CW] = noise_generator(rho_values(i), N_values(j));
            CW_est = cov(w);
            diff = CW - CW_est;
            mse_temp(k) = mean(diff(:).^2);
        end
        MSE(i,j) = mean(mse_temp);  % average over Monte Carlo runs
    end
end
MSE_dB = 10 * log10(MSE);
%% --- Plot 1: MSE (dB) vs rho for different N ---
figure;
hold on;

for i = 1:length(N_values)
    plot(rho_values, MSE_dB(:, i), '-o', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('N = %d', N_values(i)));
end

xlabel('\rho');
ylabel('Mean Squared Error (dB)');
title('MSE vs \rho for different N');
legend('Location', 'best');
grid on;
hold off;


% --- Plot 2: MSE (dB) vs L for different rho ---
figure;
hold on;

for j = 1:length(rho_values)
    plot(N_values, MSE_dB(j, :), '-s', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\rho = %.2f', rho_values(j)));
end

xlabel('N (Number of antennas)');
ylabel('Mean Squared Error (dB)');
title('MSE vs N for different \rho');
legend('Location', 'best');
grid on;
hold off;
