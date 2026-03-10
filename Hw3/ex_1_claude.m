%% =========================================================
%  HW3 – Trading & Best Combination of Stocks in Portfolio
%  AA 2025-2026
%  Data source: DataStocks/Hw3.mat
%  Each stock variable is 1326x5: [Open, High, Low, Close, Volume]
%  'Stock' = 10x4 char (ticker names)
%  'T'     = 5x6  char (time/date labels)
%
%  OLS linear predictor with sector-aware features.
%  Sectors (from HW3 spec):
%    Social    : FB, TWTR
%    Software  : MSFT, IBM, SNE
%    Technology: QCOM, INTC, AVGO
%    Online    : AMZN
%    Phones    : AAPL
% ==========================================================
clear; clc; close all;

%% ---- 0. LOAD DATA ------------------------------------------
load("DataStocks/Hw3.mat");

% Parse ticker names from 'Stock' (10x4 char array)
TICKERS  = strtrim(cellstr(Stock));
N_STOCKS = numel(TICKERS);

fprintf('Loaded %d stocks: ', N_STOCKS);
fprintf('%s  ', TICKERS{:}); fprintf('\n');

% Column layout of each 1326×5 matrix
COL_OPEN   = 1;
COL_HIGH   = 2;
COL_LOW    = 3;
COL_CLOSE  = 4;
COL_VOLUME = 5;

% Stack into [T_total × N_STOCKS] matrices
tmp     = eval(TICKERS{1});
T_total = size(tmp,1);

Open   = zeros(T_total, N_STOCKS);
High   = zeros(T_total, N_STOCKS);
Low    = zeros(T_total, N_STOCKS);
Close  = zeros(T_total, N_STOCKS);
Volume = zeros(T_total, N_STOCKS);

for k = 1:N_STOCKS
    M = eval(TICKERS{k});
    Open(:,k)   = M(:,COL_OPEN);
    High(:,k)   = M(:,COL_HIGH);
    Low(:,k)    = M(:,COL_LOW);
    Close(:,k)  = M(:,COL_CLOSE);
    Volume(:,k) = M(:,COL_VOLUME);
end
fprintf('Total trading days in .mat: %d\n', T_total);

%% ---- 1. SECTOR DEFINITION ----------------------------------
% Assign each ticker to a sector index.
% Sector mean intraday return will be added as an extra feature
% for each stock's predictor (sector momentum).
%
%  1 = Social     (FB, TWTR)
%  2 = Software   (MSFT, IBM, SNE)
%  3 = Technology (QCOM, INTC, AVGO)
%  4 = Online     (AMZN)
%  5 = Phones     (AAPL)
%  0 = Other/Unknown (ADBE, AMD, GOOGL – not in spec; treated as own sector)

SECTOR_NAMES = {'Social','Software','Technology','Online','Phones','Other'};

SECTOR_MAP = containers.Map( ...
    {'FB','TWTR',  'MSFT','IBM','SNE',  'QCOM','INTC','AVGO',  'AMZN',  'AAPL', ...
     'ADBE','AMD','GOOGL'}, ...
    {  1,    1,      2,    2,    2,       3,     3,     3,       4,       5, ...
       6,    6,    6} );

sector_idx = zeros(N_STOCKS,1);
for k = 1:N_STOCKS
    if isKey(SECTOR_MAP, TICKERS{k})
        sector_idx(k) = SECTOR_MAP(TICKERS{k});
    else
        sector_idx(k) = 6;   % Other
    end
    fprintf('  %-6s -> sector %d (%s)\n', TICKERS{k}, sector_idx(k), ...
            SECTOR_NAMES{sector_idx(k)});
end

% For convenience: list of stock indices belonging to each sector
N_SECTORS = 6;
sector_members = cell(N_SECTORS,1);
for s = 1:N_SECTORS
    sector_members{s} = find(sector_idx == s);
end

%% ---- 2. EXTRACT 2018 WINDOW --------------------------------
% 1326 days ≈ Jan-2013 to Dec-2018  (252 trading days/year).
% 2018 starts at row:  5*252+1 = 1261
% Adjust YEAR_START if your data begins at a different year.
DAYS_PER_YEAR = 252;
YEAR_START    = 2013;   % <-- adjust if needed
TARGET_YEAR   = 2018;

idx_2018_start = (TARGET_YEAR - YEAR_START) * DAYS_PER_YEAR + 1;
idx_2018_end   = min(idx_2018_start + DAYS_PER_YEAR - 1, T_total);
T_2018         = idx_2018_end - idx_2018_start + 1;

% Training = first half of 2018 (Jan–Jun), Testing = second half (Jul–Dec)
T_train = floor(T_2018 / 2);
T_test  = T_2018 - T_train;

fprintf('\n2018 rows %d–%d  (%d days) | Train=%d | Test=%d\n', ...
        idx_2018_start, idx_2018_end, T_2018, T_train, T_test);

% Include LAG rows before 2018 for warm-up
LAG    = 5;
HIST   = max(1, idx_2018_start - LAG);
OFFSET = idx_2018_start - HIST;   % number of pre-2018 warm-up rows

% Extended arrays [warm-up + 2018]
Open_e   = Open  (HIST:idx_2018_end, :);
Close_e  = Close (HIST:idx_2018_end, :);
High_e   = High  (HIST:idx_2018_end, :);
Low_e    = Low   (HIST:idx_2018_end, :);
Volume_e = Volume(HIST:idx_2018_end, :);
T_e      = size(Open_e,1);

% 2018-only slices
Open_y   = Open_e  (OFFSET+1:end, :);
Close_y  = Close_e (OFFSET+1:end, :);

%% ---- 3. DERIVED SERIES (on extended array) -----------------
% Volume normalised by 20-day rolling mean
Vol_n = zeros(size(Volume_e));
for k = 1:N_STOCKS
    for t = 1:T_e
        w = max(1,t-19):t;
        Vol_n(t,k) = Volume_e(t,k) / (mean(Volume_e(w,k)) + eps);
    end
end

% Intraday return  (close-open)/open  — proxy for Pk(t+1|t) vs Pk(t|t)
Intra_e = (Close_e - Open_e) ./ (Open_e + eps);

% Normalised high-low range  (High-Low)/Open
HL_e = (High_e - Low_e) ./ (Open_e + eps);

% Sector mean intraday return  [T_e x N_SECTORS]
% For each day t and sector s: mean of Intra_e(t, members of s)
Sector_intra = zeros(T_e, N_SECTORS);
for s = 1:N_SECTORS
    mem = sector_members{s};
    if ~isempty(mem)
        Sector_intra(:,s) = mean(Intra_e(:,mem), 2);
    end
end

%% ---- 4. OLS LINEAR PREDICTOR TRAINING ----------------------
% Feature vector x_k(t) [at extended index te]:
%
%  Group A – own-stock features (lags 1..LAG):
%    close_k(te-j), open_k(te-j), vol_norm_k(te-j),
%    intra_ret_k(te-j), HL_range_k(te-j)
%
%  Group B – current day:
%    open_k(te)  [known at decision time]
%
%  Group C – cross-stock / market (lag 1):
%    mean intraday return of ALL stocks at te-1  (market momentum)
%
%  Group D – sector features (lag 1):
%    mean intraday return of the stock's OWN sector at te-1
%    mean intraday return of each OTHER sector at te-1  (N_SECTORS-1 values)
%    → gives the predictor information about inter-sector dynamics
%
%  Bias term: 1
%
%  Target: close_k(te)  = P_k(t+1|t)

B     = cell(N_STOCKS,1);
R2_tr = zeros(N_STOCKS,1);

fprintf('\n--- OLS Training ---\n');
for k = 1:N_STOCKS
    sk = sector_idx(k);   % this stock's sector
    Xmat = [];
    yvec = [];

    for t = 1:T_train
        te = t + OFFSET;
        if te <= LAG, continue; end
        x  = feat(te, k, sk, Open_e, Close_e, Vol_n, Intra_e, HL_e, ...
                  Sector_intra, LAG, N_SECTORS);
        Xmat = [Xmat; x];
        yvec = [yvec; Close_y(t,k)];
    end

    % Pure OLS: b = (X'X)^{-1} X'y
    B{k}     = (Xmat' * Xmat) \ (Xmat' * yvec);
    yhat     = Xmat * B{k};
    ss_res   = sum((yvec - yhat).^2);
    ss_tot   = sum((yvec - mean(yvec)).^2);
    R2_tr(k) = 1 - ss_res / ss_tot;

    fprintf('  %-6s  sector=%-12s  R²_train=%.4f  #features=%d\n', ...
            TICKERS{k}, SECTOR_NAMES{sk}, R2_tr(k), numel(B{k}));
end

%% ---- 5. TRADING SIMULATION (full 2018) ---------------------
% Decision rule each day t:
%   Compute predicted intraday gain for every stock:
%     g_pred_k = ( P̂_k(t+1|t) - P_k(t|t) ) / P_k(t|t)
%   Invest ALL capital in the stock with the highest g_pred_k, IF g_pred_k > 0.
%   If all predictions are negative → hold cash (no position).

INITIAL_CAPITAL = 10000;   % €
capital         = INITIAL_CAPITAL;
cap_series      = nan(T_2018,1);
chosen          = zeros(T_2018,1);      % 0=cash, 1..N = stock index
g_pred_all      = zeros(T_2018, N_STOCKS);
g_actual_all    = zeros(T_2018, N_STOCKS);

for t = 1:T_2018
    te = t + OFFSET;
    if te <= LAG
        cap_series(t) = capital;
        continue;
    end

    % Predict closing price for each stock
    pc = zeros(N_STOCKS,1);
    for k = 1:N_STOCKS
        sk = sector_idx(k);
        x     = feat(te, k, sk, Open_e, Close_e, Vol_n, Intra_e, HL_e, ...
                     Sector_intra, LAG, N_SECTORS);
        pc(k) = x * B{k};
    end

    % Predicted and actual intraday % gain
    gp = (pc            - Open_y(t,:)') ./ (Open_y(t,:)' + eps);
    ga = (Close_y(t,:)  - Open_y(t,:))  ./ (Open_y(t,:)  + eps);
    g_pred_all(t,:)   = gp';
    g_actual_all(t,:) = ga;

    % Invest all capital in best predicted stock (if gain > 0)
    [best_g, best_k] = max(gp);
    if best_g > 0
        capital   = (capital / Open_y(t,best_k)) * Close_y(t,best_k);
        chosen(t) = best_k;
    else
        chosen(t) = 0;   % hold cash
    end
    cap_series(t) = capital;
end

%% ---- 6. RESULTS --------------------------------------------
fprintf('\n========== RESULTS ==========\n');
fprintf('Initial capital  : €%10.2f\n', INITIAL_CAPITAL);
fprintf('Final   capital  : €%10.2f\n', capital);
fprintf('Total return     : %+.2f%%\n',  (capital/INITIAL_CAPITAL-1)*100);

fprintf('\nBuy-and-hold 2018 benchmarks:\n');
for k = 1:N_STOCKS
    bh = INITIAL_CAPITAL * Close_y(end,k) / Open_y(1,k);
    fprintf('  %-6s  €%9.2f  (%+.1f%%)\n', TICKERS{k}, bh, (bh/INITIAL_CAPITAL-1)*100);
end

dC  = diff(cap_series);
dC(isnan(dC)) = [];
pos = sum(dC > 0);
neg = sum(dC <= 0);
fprintf('\nMean daily ΔC  : €%+.2f\n',   mean(dC));
fprintf('Std  daily ΔC  : €%.2f\n',       std(dC));
fprintf('Sharpe ratio   : %.3f\n',         mean(dC)/std(dC)*sqrt(252));
fprintf('Win rate       : %d/%d (%.1f%%)\n', pos, pos+neg, 100*pos/(pos+neg));

fprintf('\nDays invested per stock:\n');
for k = 1:N_STOCKS
    fprintf('  %-6s  %3d days\n', TICKERS{k}, sum(chosen==k));
end
fprintf('  %-6s  %3d days\n', 'Cash', sum(chosen==0));

fprintf('\nDirectional accuracy on TEST period:\n');
for k = 1:N_STOCKS
    sk = sector_idx(k);
    cor = 0; tot = 0;
    for t = T_train+1 : T_2018
        te = t + OFFSET;
        if te <= LAG, continue; end
        x   = feat(te, k, sk, Open_e, Close_e, Vol_n, Intra_e, HL_e, ...
                   Sector_intra, LAG, N_SECTORS);
        pred_up   = (x * B{k}) > Open_y(t,k);
        actual_up = Close_y(t,k) > Open_y(t,k);
        cor = cor + (pred_up == actual_up);
        tot = tot + 1;
    end
    fprintf('  %-6s  %.1f%%  (%d/%d)\n', TICKERS{k}, 100*cor/tot, cor, tot);
end

%% ---- 7. FIGURES --------------------------------------------
days = (1:T_2018)';

% -- Fig 1: Capital evolution
figure('Name','Capital Evolution','Position',[50 50 1100 420]);
plot(days, cap_series, 'b-', 'LineWidth', 1.8); hold on;
xline(T_train, 'r--', 'LineWidth', 1.5, 'Label', 'Train | Test');
yline(INITIAL_CAPITAL, 'k:', 'LineWidth', 1.2, 'Label', '€10,000');
xlabel('Trading day (2018)'); ylabel('Portfolio value (€)');
title('Portfolio Capital – OLS Linear Predictor, Single-Stock Strategy');
legend('Capital','Train/Test split','Initial capital','Location','northwest');
grid on; set(gca,'FontSize',11);

% -- Fig 2: Normalised stock prices (colour-coded by sector)
sector_colors = lines(N_SECTORS);
figure('Name','Normalised Prices','Position',[50 500 1100 400]);
p_norm = Close_y ./ Close_y(1,:);
hold on;
for k = 1:N_STOCKS
    plot(days, p_norm(:,k), 'Color', sector_colors(sector_idx(k),:), 'LineWidth', 1.3);
end
xlabel('Trading day (2018)'); ylabel('Price / Price_{day 1}');
title('2018 Normalised Stock Prices  (colour = sector)');
% Custom legend: one entry per sector
h = gobjects(N_SECTORS,1);
for s = 1:N_SECTORS
    h(s) = plot(nan, nan, '-', 'Color', sector_colors(s,:), 'LineWidth', 2);
end
legend(h, SECTOR_NAMES, 'Location','best');
grid on; set(gca,'FontSize',11);

% -- Fig 3: Stock selection timeline
figure('Name','Stock Selection','Position',[50 50 1100 280]);
scatter(days, chosen, 25, chosen, 'filled');
yticks(0:N_STOCKS);
yticklabels(['Cash'; TICKERS(:)]);
xlabel('Trading day (2018)');
title('Stock Selected Each Day  (0 = Cash)');
colormap(lines(N_STOCKS+1)); grid on; set(gca,'FontSize',10);

% -- Fig 4: Predicted gain heatmap
figure('Name','Predicted Gain Heatmap','Position',[50 50 1100 370]);
imagesc(days, 1:N_STOCKS, g_pred_all');
cmap = redblue(64); colormap(cmap); colorbar;
clim([-0.03 0.03]);
set(gca,'YTick',1:N_STOCKS,'YTickLabel',TICKERS,'FontSize',10);
xlabel('Trading day (2018)');
title('Predicted Intraday Gain per Stock  (red > 0, blue < 0)');

% -- Fig 5: Daily ΔC histogram
figure('Name','Daily P&L','Position',[200 200 620 370]);
histogram(dC, 40, 'FaceColor',[0.2 0.5 0.9], 'EdgeColor','none'); hold on;
xline(0,     'r-',  'LineWidth', 2);
xline(mean(dC), 'g--','LineWidth', 1.5, 'Label', sprintf('Mean=€%.1f',mean(dC)));
xlabel('Daily ΔC (€)'); ylabel('Count');
title(sprintf('Daily P&L Distribution  (mean=€%.1f,  σ=€%.1f)', mean(dC), std(dC)));
grid on; set(gca,'FontSize',11);

% -- Fig 6: OLS R² per stock (colour-coded by sector)
figure('Name','R² per Stock','Position',[200 200 720 370]);
bar_colors = sector_colors(sector_idx,:);
b = bar(R2_tr, 'FaceColor','flat');
b.CData = bar_colors;
set(gca,'XTick',1:N_STOCKS,'XTickLabel',TICKERS,'FontSize',11);
ylabel('R²  (training set)'); title('OLS Predictor Training R²  (colour = sector)');
ylim([0 1]); grid on;
h2 = gobjects(N_SECTORS,1);
for s = 1:N_SECTORS
    h2(s) = patch(nan,nan, sector_colors(s,:));
end
legend(h2, SECTOR_NAMES, 'Location','northeast');

% -- Fig 7: Sector mean intraday return (2018 window)
Sector_intra_y = Sector_intra(OFFSET+1:end, :);
figure('Name','Sector Returns','Position',[200 200 1100 380]);
plot(days, cumsum(Sector_intra_y), 'LineWidth', 1.5);
legend(SECTOR_NAMES, 'Location','best');
xlabel('Trading day (2018)'); ylabel('Cumulative intraday return');
title('Cumulative Sector Returns (2018)'); grid on; set(gca,'FontSize',11);

fprintf('\n--- Simulation complete ---\n');

%% =========================================================
%  LOCAL FUNCTIONS
% ==========================================================
function x = feat(te, k, sk, Open_e, Close_e, Vol_n, Intra_e, HL_e, ...
                  Sector_intra, L, N_SECTORS)
%FEAT  Build OLS feature row for stock k (sector sk) at extended index te.
%
%  Features:
%   [1]           bias
%   [1]           current open  P_k(te)                     (Group B)
%   [5*L]         lagged own-stock: close, open, vol_norm,   (Group A)
%                   intra_ret, HL_range  for lags 1..L
%   [1]           market momentum: mean intra of ALL stocks at te-1  (Group C)
%   [1]           own-sector momentum: mean intra of sector sk at te-1 (Group D)
%   [N_SECTORS-1] other-sector momenta at te-1               (Group D)

    % Group B: current open
    f = [1, Open_e(te,k)];

    % Group A: lagged own-stock features (lags 1..L)
    for lag = 1:L
        c  = Close_e(te-lag, k);
        o  = Open_e(te-lag, k);
        v  = Vol_n(te-lag, k);
        ir = Intra_e(te-lag, k);
        hl = HL_e(te-lag, k);
        f  = [f, c, o, v, ir, hl];
    end

    % Group C: market-wide momentum (mean intraday return, lag 1)
    f = [f, mean(Intra_e(te-1,:))];

    % Group D: sector momenta (lag 1)
    % own sector return, then each other sector's return
    own_ret   = Sector_intra(te-1, sk);
    all_sec   = Sector_intra(te-1, :);
    other_idx = [1:sk-1, sk+1:N_SECTORS];
    other_ret = all_sec(other_idx);
    f = [f, own_ret, other_ret];

    x = f;
end

function cmap = redblue(n)
%REDBLUE  Diverging blue-white-red colormap with n levels
    h = floor(n/2);
    r = [linspace(0,1,h),   ones(1,n-h)];
    g = [linspace(0,1,h),   linspace(1,0,n-h)];
    b = [ones(1,h),          linspace(1,0,n-h)];
    cmap = [r', g', b'];
end