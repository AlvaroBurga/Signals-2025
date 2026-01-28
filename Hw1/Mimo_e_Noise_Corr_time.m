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

%% parameters

betta = 0.5;
N = 4;
M = 10^3;
L = 5;

alfa = 0.1;
SNR = 100;
rho = 0.999;

%% Real values initialization
function [h, h_est, CRB] = MIMO_Indentification(SNR, rho, alfa)
    variance_w = 1/SNR;

    P = M + L - 1;
    [w,Cw_antenna] = noise_generator(rho, N, variance_w, P);
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
    
    CRB = inv(X' * X) * variance_w; %Equal in all y, it only depends on X.
end

%% Error estimation
rho_values = [0, 0.5, 0.9];
alfa_values = [0.1, 0.9];
SNR_dB_values = [-10, -5, 0, 5, 10, 15, 20, 30, 40];  % SNR in dB
SNR_values = 10.^(SNR_dB_values / 10);

diff = h - h_est;
MSE = mean(diff(:).^2)
MSE_dB = 10*log10(MSE)