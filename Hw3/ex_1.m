%% H3 Alvaro Burga
clc;
clear;
close all;

%% Load data
load("DataStocks/Hw3.mat");

open_prices   = table();
close_prices  = table();
total_volumes = table();

for i = 1:length(Stock)
    ticker = strtrim(Stock(i,:));
    ticker_data = eval(ticker);
    open_prices.(ticker)   = ticker_data(2:end, 1);
    close_prices.(ticker)  = ticker_data(2:end, 4);
    total_volumes.(ticker) = ticker_data(2:end, 5);
end

T_days   = height(open_prices);
n_stocks = width(open_prices);
tickers  = open_prices.Properties.VariableNames;

P_o = table2array(open_prices);
P_c = table2array(close_prices);

% Jan-Jun: ~124 trading days | Jul-Dec: rest
N_train = 124;
idx_tr  = 1:N_train;
idx_te  = (N_train+1):T_days;

%% Train one linear predictor per stock
% x_k(t) = [1, open(t), close(t-1), close(t-2), close(t-3), close(t-4)]
L = 5;
B = zeros(L+1, n_stocks);

for k = 1:n_stocks
    valid = idx_tr(idx_tr > L);
    N_tr  = length(valid);
    Xtr   = zeros(N_tr, L+1);
    ytr   = zeros(N_tr, 1);
    for ii = 1:N_tr
        t = valid(ii);
        Xtr(ii,:) = [1, P_o(t,k), P_c(t-1,k), P_c(t-2,k), P_c(t-3,k), P_c(t-4,k)];
        ytr(ii)   = P_c(t,k);
    end
    B(:,k) = (Xtr'*Xtr) \ (Xtr'*ytr);
end

%% Simulate trading on test period
C0 = 10000;
C  = zeros(length(idx_te)+1, 1);
C(1) = C0;
chosen      = zeros(length(idx_te), 1);
pred_gain_m = zeros(length(idx_te), n_stocks);
act_gain_m  = zeros(length(idx_te), n_stocks);

for ii = 1:length(idx_te)
    t = idx_te(ii);

    P_hat = zeros(n_stocks, 1);
    for k = 1:n_stocks
        x_t = [1, P_o(t,k), P_c(t-1,k), P_c(t-2,k), P_c(t-3,k), P_c(t-4,k)]';
        P_hat(k) = B(:,k)' * x_t;
    end

    pred_gain_m(ii,:) = (P_hat' - P_o(t,:)) ./ P_o(t,:);
    act_gain_m(ii,:)  = (P_c(t,:)  - P_o(t,:)) ./ P_o(t,:);

    [best_g, best_k] = max(pred_gain_m(ii,:));
    chosen(ii) = best_k;

    if best_g > 0
        N_sh = C(ii) / P_o(t, best_k);
        C(ii+1) = N_sh * P_c(t, best_k);
    else
        C(ii+1) = C(ii);
    end
end

fprintf('Initial capital:  %.2f EUR\n', C0);
fprintf('Final capital:    %.2f EUR\n', C(end));
fprintf('Total gain:       %.2f%%\n', 100*(C(end)-C0)/C0);

%% --- Plot 1: Capital evolution ---
figure;
plot(1:length(C), C, 'LineWidth', 1.5);
xlabel('Trading day (Jul-Dec 2018)');
ylabel('Capital (EUR)');
title('Capital evolution - linear predictor strategy');
grid on;

%% --- Plot 2: Stock selection frequency ---
figure;
histogram(chosen(chosen > 0), 0.5:1:n_stocks+0.5);
xticks(1:n_stocks);
xticklabels(tickers);
xlabel('Stock');
ylabel('Days selected');
title('Days each stock was chosen as best');
grid on;

%% --- Plot 3: Directional accuracy per stock ---
correct = zeros(1, n_stocks);
for ii = 1:length(idx_te)
    actual_dir = act_gain_m(ii,:) > 0;
    pred_dir   = pred_gain_m(ii,:) > 0;
    correct    = correct + (actual_dir == pred_dir);
end
accuracy = correct / length(idx_te) * 100;

figure;
bar(accuracy);
xticks(1:n_stocks);
xticklabels(tickers);
ylabel('Accuracy (%)');
title('Directional accuracy of linear predictor per stock');
ylim([0 100]);
grid on;
