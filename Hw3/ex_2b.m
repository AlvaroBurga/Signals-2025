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

%% Train logistic models (first semester)
L = 5;
sigmoid = @(z) 1 ./ (1 + exp(-z));

B_log   = zeros(L+2, n_stocks);
mu_all  = zeros(n_stocks, L+1);
sig_all = zeros(n_stocks, L+1);

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
    mu_x  = mean(Xtr(:,2:end));
    sig_x = std(Xtr(:,2:end));
    Xtr_n = [ones(N_tr,1), (Xtr(:,2:end) - mu_x) ./ sig_x];

    b = zeros(L+2, 1);
    for it = 1:2000
        hh   = sigmoid(Xtr_n * b);
        grad = Xtr_n' * (hh - ytr);
        b    = b - (0.1 / N_tr) * grad;
    end
    B_log(:,k)   = b;
    mu_all(k,:)  = mu_x;
    sig_all(k,:) = sig_x;
end

%% Simulate trading - 3 strategies
% S1: winner-takes-all  (all capital on stock with max h, only if h > 0.5)
% S2: proportional      (capital split proportional to h_k for h_k > 0.5)
% S3: high-confidence   (capital split proportional to h_k for h_k > 0.6)

C0 = 10000;
N_te = length(idx_te);
C1 = zeros(N_te+1, 1);  C1(1) = C0;
C2 = zeros(N_te+1, 1);  C2(1) = C0;
C3 = zeros(N_te+1, 1);  C3(1) = C0;

for ii = 1:N_te
    t = idx_te(ii);

    % Compute h_k for all stocks
    h = zeros(n_stocks, 1);
    for k = 1:n_stocks
        x_t = [1, (P_o(t,k)       - mu_all(k,1))/sig_all(k,1), ...
                   (P_c(t-1,k)    - mu_all(k,2))/sig_all(k,2), ...
                   (P_c(t-2,k)    - mu_all(k,3))/sig_all(k,3), ...
                   (P_c(t-3,k)    - mu_all(k,4))/sig_all(k,4), ...
                   (P_c(t-4,k)    - mu_all(k,5))/sig_all(k,5), ...
                   (Vol(t-1,k)/1e6 - mu_all(k,6))/sig_all(k,6)]';
        h(k) = sigmoid(B_log(:,k)' * x_t);
    end

    daily_ret = P_c(t,:)' ./ P_o(t,:)';  % actual return ratio each stock

    % S1: winner-takes-all
    [best_h, best_k] = max(h);
    if best_h > 0.5
        C1(ii+1) = C1(ii) * daily_ret(best_k);
    else
        C1(ii+1) = C1(ii);
    end

    % S2: proportional to h (only h > 0.5)
    mask2 = h > 0.5;
    if any(mask2)
        w2 = h .* mask2;
        w2 = w2 / sum(w2);
        C2(ii+1) = C2(ii) * (w2' * daily_ret);
    else
        C2(ii+1) = C2(ii);
    end

    % S3: high-confidence (h > 0.6)
    mask3 = h > 0.6;
    if any(mask3)
        w3 = h .* mask3;
        w3 = w3 / sum(w3);
        C3(ii+1) = C3(ii) * (w3' * daily_ret);
    else
        C3(ii+1) = C3(ii);
    end
end

fprintf('Strategy             Final capital   Gain\n');
fprintf('S1 winner-takes-all  %10.2f EUR  %+.1f%%\n', C1(end), 100*(C1(end)-C0)/C0);
fprintf('S2 proportional      %10.2f EUR  %+.1f%%\n', C2(end), 100*(C2(end)-C0)/C0);
fprintf('S3 high-confidence   %10.2f EUR  %+.1f%%\n', C3(end), 100*(C3(end)-C0)/C0);

%% --- Plot 1: Capital evolution - all strategies ---
t_ax = 0:N_te;
figure;
hold on;
plot(t_ax, C1, '-',  'LineWidth', 1.5, 'DisplayName', 'S1: winner-takes-all');
plot(t_ax, C2, '--', 'LineWidth', 1.5, 'DisplayName', 'S2: proportional (h>0.5)');
plot(t_ax, C3, ':',  'LineWidth', 2,   'DisplayName', 'S3: high-confidence (h>0.6)');
yline(C0, 'k--', 'LineWidth', 1, 'DisplayName', 'initial capital');
xlabel('Trading day (Jul-Dec 2018)');
ylabel('Capital (EUR)');
title('Portfolio capital - logistic regression strategies');
legend('Location', 'best');
grid on;
hold off;

%% --- Plot 2: Daily h values heatmap over test period ---
H_test = zeros(N_te, n_stocks);
for ii = 1:N_te
    t = idx_te(ii);
    for k = 1:n_stocks
        x_t = [1, (P_o(t,k)       - mu_all(k,1))/sig_all(k,1), ...
                   (P_c(t-1,k)    - mu_all(k,2))/sig_all(k,2), ...
                   (P_c(t-2,k)    - mu_all(k,3))/sig_all(k,3), ...
                   (P_c(t-3,k)    - mu_all(k,4))/sig_all(k,4), ...
                   (P_c(t-4,k)    - mu_all(k,5))/sig_all(k,5), ...
                   (Vol(t-1,k)/1e6 - mu_all(k,6))/sig_all(k,6)]';
        H_test(ii,k) = sigmoid(B_log(:,k)' * x_t);
    end
end

figure;
imagesc(H_test');
colorbar;
clim([0 1]);
yticks(1:n_stocks);
yticklabels(tickers);
xlabel('Trading day (test)');
title('Sigmoid output h_k(t) - test period');
colormap(jet);

%% --- Plot 3: Final capital comparison ---
figure;
bar([C1(end), C2(end), C3(end), C0]);
xticklabels({'S1 winner', 'S2 proportional', 'S3 high-conf', 'Initial'});
ylabel('Capital (EUR)');
title('Final capital by strategy (Jul-Dec 2018)');
grid on;
