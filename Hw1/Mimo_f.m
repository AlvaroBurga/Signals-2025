%% H1 Alvaro Burga
clc;
clear;
close all;

%% Function - noise
function [w, CW] = noise_generator(rho, N,variance_w, P)

    CW = zeros(N,N);
    for j = 1:N
        for i = 1:N
            CW(i,j) = variance_w*rho^(abs(i-j));
        end
    end
    [V,D] = eig(CW);
    w = zeros(P,N);
    for k = 1:P
        z = randn(N,1);
        w(k,:) = V* sqrt(D) * z;
    end
end


%% Identification
function [h_est, y, x] = MIMO_Indentification(SNR, rho, alfa)
    variance_w = 1/SNR;
    betta = 0.5;
    N = 4;
    M = 10^3;
    L = 5;

    P = M + L - 1;
    [w,~] = noise_generator(rho, N, variance_w, P);
    x = randn(M, N);
    h = zeros(N, N,L);
    
    %Generate matrix h
    for i = 1:N
        for j = 1:N
            for k = 1:L
                h(i,j,k) = alfa^(abs(i-j)) * betta^k;
            end
        end
    end
    
    %Generate y response
    y = zeros (P,N);
    for i = 1:N
        y_clean = zeros(P,1);
        for j = 1:N
            hij = reshape(h(i, j, :), 5, 1);
            H = convmtx(hij,M);
            y_clean = y_clean + H*x(:,j) ;
        end
        y(:,i) = y_clean + w(:,i);
    end
    
    % Identification
    h_est = zeros(N,N,L);
    X = zeros (P,N*L);
    
    for j = 1:N %For each antenna
        f=(j-1)*L;
        X(:,1+f:L+f) = convmtx(x(:,j),L);
    end
    
    for i = 1:N %For each y
    
        h_est_aux = inv(X'*X)*X'*y(:,i);
        for j = 1:N %For each antenna
            f=(j-1)*L;
            h_est(i,j,:) = h_est_aux(1+f:L+f);
        end
    end
   
end
%% Deconvolution  (RAM Killer)
function [x_real, x_est, CRB] = MIMO_Deconvolution(SNR, rho, alfa)

    [h_est,y, x_real] = MIMO_Indentification(SNR, rho, alfa);
    variance_w = 1/SNR;
    N = 4;
    M = 10^3;
    L = 5;
    P = M + L - 1;
    y_hat = y(:);
    H_hat = zeros(N*P, N*M);

    for i = 1:N
        for j = 1:N
            hij = reshape(h_est(i, j, :), 5, 1);
            Hij = convmtx(hij,M);
            H_hat((i-1)*P+1:i*P, (j-1)*M+1:j*M) = Hij;
        end
    end
    I = eye(4*P);
    x_hat = H_hat' * inv(H_hat *H_hat' + variance_w*I)*y_hat;
    CRB = inv(H_hat' *H_hat/variance_w + eye(4*M));
    x_est = zeros(M,N);
    for j = 1:N
        x_est(:,j) = x_hat((j-1)*M+1:j*M);
    end
end

%% Error estimation
%rho_values = [0, 0.5, 0.9];
rho_values = [0.5];
alfa_values = [0.1, 0.9];
SNR_dB_values = [-40,-10, 0, 10, 20, 30, 40, 80];  % SNR in dB
SNR_values = 10.^(SNR_dB_values / 10);

% Preallocate MSE matrix
MC = 5;  % number of Monte Carlo runs
MSE = zeros(length(SNR_values), length(alfa_values), length(rho_values));
CRB_values = zeros(length(SNR_values), length(alfa_values), length(rho_values));
cont = 1;
for i = 1:length(SNR_values)
    for j = 1:length(alfa_values)
        for k = 1:length(rho_values) 
            mse_temp = zeros(1, MC);  % store MSE for each run
            CRB_temp = zeros(1, MC);
            for m = 1:MC
                [x_real, x_est,CRB_matrix] = MIMO_Deconvolution(SNR_values(i),rho_values(k), alfa_values(j));
                diff = x_real - x_est;
                mse_temp(m) = mean(diff(:).^2);
                CRB_temp(m) = mean(diag(CRB_matrix));
                
            end
            progress = cont*100/(length(SNR_values) * length(alfa_values) * length(rho_values));
            cont = cont + 1;
            fprintf('Progress: %.2f%%\n', progress);
            MSE(i,j,k) = mean(mse_temp);  % average over Monte Carlo runs
            CRB_values(i,j,k) = mean(CRB_temp);
        end
    end

end
MSE_dB = 10 * log10(MSE);
CRB_dB = 10 * log10(CRB_values);

%% --- Plot 1: MSE (dB) vs SNR (dB) for different alfa ---
figure;
hold on;
for i = 1:length(alfa_values)
    plot(SNR_dB_values, MSE_dB(:,i,1), '-o', 'LineWidth', 1.5, ...
        'DisplayName', ['alfa = ' num2str(alfa_values(i))]);
end
xlabel('SNR (dB)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR for different alfa (fixed \rho)');
legend('Location', 'northeast');
grid on;
hold off;

% --- Plot 2: MSE (dB) vs SNR (dB) for different rho ---
figure;
hold on;
for j = 1:length(rho_values)
    plot(SNR_dB_values, MSE_dB(:,1,j), '-o', 'LineWidth', 1.5, ...
        'DisplayName', ['\rho = ' num2str(rho_values(j))]);
end
xlabel('SNR (dB)');
ylabel('Mean Squared Error (dB)');
title('MSE vs SNR for different \rho (fixed alfa)');
legend('Location', 'northeast');
grid on;
hold off;

% --- Plot 3: MSE (dB) vs SNR (dB) vs CRB(dB)
figure;
hold on;
plot(SNR_dB_values, MSE_dB(:,1,1), '-o', 'LineWidth', 1.5, ...
    'DisplayName', ['MSE']);
plot(SNR_dB_values, CRB_dB(:,1,1), '-o', 'LineWidth', 1.5, ...
    'DisplayName', ['CRB']);
xlabel('SNR (dB)');
ylabel('Value (dB)');
title('MSE and CRB vs SNR');
legend('Location', 'northeast');
grid on;
hold off;
