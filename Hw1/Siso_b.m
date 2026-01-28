%% Lab 2, montecarlo
clc;
clear;
close all;

%% SISO b function
function [h, y, h_real] = h_identification_with_correlated_noise(L , SNR)
    N = 200;
    M = L + N - 1;
    variance_w = 1 / SNR;

    b = 1;            % numerator
    a = [1 -0.9];     % denominator
    delta = [1; zeros(L-1,1)];  % unit impulse
    h_real = filter(b, a, delta);
    
    % We create the correlation noise variables
    [i, j] = meshgrid(1:M, 1:M); 
    Cw = variance_w *(-1/9).^abs(i - j); 
    CwInv = inv(Cw);

    % We use this artifcat to find and equivalence noise from the
    % Correlation
    [V,D] = eig(Cw);
    z = randn(M,1);
    w = V* sqrt(D) * z;
    
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
    
    h = inv(X'*CwInv* X)* X' * CwInv * y;

    
end

%% SISO b results

SNR_dB_values = [-10, -5, 0, 5, 10, 15, 20, 30, 40];  % SNR in dB
SNR_values = 10.^(SNR_dB_values / 10);
L_values = [2,5,10,20,50,100];
MC = 100;  % number of Monte Carlo runs


% Preallocate MSE matrix
MSE = zeros(length(L_values), length(SNR_values));

for i = 1:length(L_values)
    for j = 1:length(SNR_values)
        mse_temp = zeros(1, MC);  % store MSE for each run
        for k = 1:MC
            [h, y, h_real] = h_identification_with_correlated_noise(L_values(i), SNR_values(j));
            mse_temp(k) = mean((h_real - h).^2);
        end
        MSE(i,j) = mean(mse_temp);  % average over Monte Carlo runs
    end
end
MSE_dB = 10 * log10(MSE);

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