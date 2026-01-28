%% Lab 2, montecarlo
clc;
clear;
close all;

%% 1. System identification and channel Estimation

h = [1 2 4 -1 2].';      %Channel impulse response
u = randn(5,1);        %Training samples

U = convmtx(u, length(u));        %Recordar para el examen. Crea la matriz conv
snr = 10;              %In db
N = length(h) + length(u) -1;

Mc = 1000;               %number of simulations in montecarlo experiment

s_w = sum(abs(h).^2) / db2pow(snr);
h_est = zeros(Mc, length(h));
for m = 1 : Mc
    x = U*h + sqrt(s_w) * randn(N,1);
    h_est(m,:) = inv(U.' * U) * U.' * x; %Recordar channel stimation
end

C_crb = s_w * inv(U.' * U) ; % El limite estadistico de la varianza del error
errbound = trace(C_crb);
mean_h = mean(h_est, 1);
% Calculate the mean squared error between estimated and true channel
mse = mean(sum((h - mean_h.').^2));

%% 1. Channel Estimation
clc;
clear;
close all;

SNR_DB_vec = -10 : 10;
Nt = length(SNR_DB_vec);

N=3;
h =  randn(N,1); 

M=3;
u =  randn(M,1);
U = convmtx(u, N);

delta = 0.01;
s_w = 1/db2pow(min(SNR_DB_vec));
L = ceil(4*s_w / delta);

for n = 1:Nt
    s_w = 1/db2pow(min(SNR_DB_vec));
    for l = 1 : L
        x = U*h + sqrt(s_w) * randn(N+M-1,1);
        h_est(l,:) = inv(U.' * U) * U.' * x; %Recordar channel stimation
    end
end