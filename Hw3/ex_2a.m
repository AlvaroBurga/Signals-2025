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
Vol = table2array(total_volumes);

N_train = 124;
idx_tr  = 1:N_train;
idx_te  = (N_train+1):T_days;

%% Logistic regression - one model per stock
% x_k(t) = [1, open(t), close(t-1), ..., close(t-4)]
L = 5;

sigmoid = @(z) 1 ./ (1 + exp(-z));

% features: [1, open(t), close(t-1..t-4), vol(t-1)/1e6]  -> L+2 total
B_log   = zeros(L+2, n_stocks);
mu_all  = zeros(n_stocks, L+1);
sig_all = zeros(n_stocks, L+1);

acc_train = zeros(1, n_stocks);
acc_test  = zeros(1, n_stocks);

h_train_all = cell(n_stocks, 1);
h_test_all  = cell(n_stocks, 1);
y_train_all = cell(n_stocks, 1);
y_test_all  = cell(n_stocks, 1);

for k = 1:n_stocks
    valid_tr = idx_tr(idx_tr > L);
    N_tr = length(valid_tr);
    Xtr  = zeros(N_tr, L+2);
    ytr  = zeros(N_tr, 1);
    for ii = 1:N_tr
        t = valid_tr(ii);
        Xtr(ii,:) = [1, P_o(t,k), P_c(t-1,k), P_c(t-2,k), P_c(t-3,k), P_c(t-4,k), Vol(t-1,k)/1e6];
        ytr(ii)   = P_c(t,k) > P_o(t,k);
    end

    % Normalize features for numerical stability
    mu_x  = mean(Xtr(:,2:end));
    sig_x = std(Xtr(:,2:end));
    Xtr_n = [ones(N_tr,1), (Xtr(:,2:end) - mu_x) ./ sig_x];

    % Gradient descent
    b     = zeros(L+2, 1);
    alpha = 0.1;
    for it = 1:2000
        h    = sigmoid(Xtr_n * b);
        grad = Xtr_n' * (h - ytr);
        b    = b - (alpha / N_tr) * grad;
    end

    B_log(:,k)    = b;
    mu_all(k,:)   = mu_x;
    sig_all(k,:)  = sig_x;

    h_tr = sigmoid(Xtr_n * b);
    acc_train(k)    = mean((h_tr > 0.5) == ytr) * 100;
    h_train_all{k}  = h_tr;
    y_train_all{k}  = ytr;

    % Test set
    valid_te = idx_te(idx_te > L);
    N_te = length(valid_te);
    Xte  = zeros(N_te, L+2);
    yte  = zeros(N_te, 1);
    for ii = 1:N_te
        t = valid_te(ii);
        Xte(ii,:) = [1, P_o(t,k), P_c(t-1,k), P_c(t-2,k), P_c(t-3,k), P_c(t-4,k), Vol(t-1,k)/1e6];
        yte(ii)   = P_c(t,k) > P_o(t,k);
    end
    Xte_n = [ones(N_te,1), (Xte(:,2:end) - mu_x) ./ sig_x];
    h_te = sigmoid(Xte_n * b);
    acc_test(k)    = mean((h_te > 0.5) == yte) * 100;
    h_test_all{k}  = h_te;
    y_test_all{k}  = yte;
end

fprintf('%-6s  Train Acc  Test Acc\n', 'Stock');
for k = 1:n_stocks
    fprintf('%-6s  %6.1f%%    %6.1f%%\n', tickers{k}, acc_train(k), acc_test(k));
end

%% --- Plot 1: Train vs test accuracy ---
figure;
bar([acc_train; acc_test]');
xticks(1:n_stocks);
xticklabels(tickers);
ylabel('Accuracy (%)');
title('Logistic regression accuracy per stock');
legend({'Train', 'Test'}, 'Location', 'southeast');
ylim([0 100]);
grid on;

%% --- Plot 2: Sigmoid output distribution (training) ---
figure;
for k = 1:n_stocks
    subplot(5,2,k);
    hold on;
    histogram(h_train_all{k}(y_train_all{k}==1), 15, 'FaceColor', [0.2 0.6 0.2]);
    histogram(h_train_all{k}(y_train_all{k}==0), 15, 'FaceColor', [0.8 0.2 0.2]);
    title(tickers{k});
    xlabel('h(x)');
    xlim([0 1]);
    hold off;
end
sgtitle('Sigmoid output - training (green=up, red=down)');

%% --- Plot 3: Sigmoid output over time - test period (first stock) ---
figure;
hold on;
plot(h_test_all{1}, 'LineWidth', 1.2, 'DisplayName', 'h(x)');
plot(y_test_all{1}, '.', 'MarkerSize', 8, 'DisplayName', 'actual label');
yline(0.5, '--', 'threshold');
xlabel('Trading day (test)');
ylabel('Probability of up');
title(['Logistic output vs actual label - ' tickers{1}]);
legend('Location', 'best');
grid on;
hold off;
