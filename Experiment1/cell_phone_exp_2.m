%% SLP 2025/2026 - Phone experiment
clc; clear; close all;
% 
% % Loading of the time and signals
load("Experiment_table.mat")

%% Initizalize

%Time and period initialization
t_all = {x1(:,1) x2(:,1) x3(:,1) x4(:,1) x5(:,1) x6(:,1)};
sig_all = {x1(:,4) x2(:,4) x3(:,4) x4(:,4) x5(:,4) x6(:,4)};

%Period and sampling
Ts = zeros(1, numel(t_all));
Fs = zeros(1, numel(t_all));
for i = 1:numel(t_all)
    ti = t_all{i};
    Ts(i) = mean(diff(ti));
    Fs(i) = 1/Ts(i);
end

%Plot in time
figure(1);
for i = 1:numel(t_all)
    plot(t_all{i},sig_all{i}); hold on;
end
legend('show');
title('Signals before synchronization');
xlabel('t[s]');

%% Syncing signals 
Fs_ref = max(Fs);

sig_resampled = cell(1,numel(sig_all)); %Cell is like a list
for i = 1:numel(sig_all)
    [P,Q] = rat(Fs_ref/Fs(i), 1e-6); %Rat is used to get the efficient int values to use resample.
    
    % Resample de la señal
    sig_resampled{i} = resample(sig_all{i}, P, Q); %Change the frecuency Q to P
end

len_all = cellfun(@length, sig_resampled);
min_len = min(len_all); 

figure(3);
for i = 1:numel(t_all)
    sig_resampled{i} = sig_resampled{i}(1:min_len)
    plot(sig_resampled{i}); hold on;
end
legend('show');
title('Signal resampled')
xlabel('Sample');
%% Signal segmentation
signals_a = cell(1, numel(sig_resampled));
signals_b = cell(1, numel(sig_resampled));

figure(4);
for i = 1:numel(sig_resampled)
    signals_a{i} = sig_resampled{i}(501:1500);
    signals_b{i} = sig_resampled{i}(1501:2500);
    plot(signals_a{i}); hold on;
end
title('Signals from experiment A')
legend('show');

%% Obtaining the sync delay between the phones
% We know that the 3-4, 2-5 and 1-6 should recieve the signal at the same
% time in the experiment A. However, due to the sync between phones, they
% should experience a time delay. Therefore, to get the time delay between
% the pairs we should take the number of samples of diference. Let's take
% in consideration that the signal 3, 2 and 1 should appear before their
% pairs in the part b.

pairs = [3 4; 2 5; 1 6];
tau_a = zeros(1, 3);
aux_a = cell(1,6);
for i = 1: 3
    figure(1); subplot(size(pairs,1), 1, i);
    plot(signals_a{pairs(i,1)}); hold on;  
    plot(signals_a{pairs(i,2)}); 
    legend(sprintf('Signal %d', pairs(i,1)), sprintf('Signal %d', pairs(i,2))); 
    xlabel('Samples'); title ('Comparing the signal pairs')

    [acor, lag] = xcorr(signals_a{pairs(i,1)},signals_a{pairs(i,2)});
    figure(2); subplot(size(pairs,1), 1, i); plot(lag,acor); 
    title ('Autocorrelation of signals');
    [~, I] = max(abs(acor));
    tau_a(i) = -lag(I);
    aux_a{pairs(i,1)} =circshift(signals_a{pairs(i,1)}, tau_a(i)); %Use this function to displace the signal tau_a
    
end

figure(3);
plot([aux_a{pairs(i,1)},2+signals_a{pairs(i,2)}],'-')
title ('Validating displacement');


%% Obtaining the propagation velocity

tau_b = zeros(1,3);
travel_space = zeros(1,3);
travel_time = zeros(1,3);
velocity = zeros(1,3);
for i = 1: 3
    x =circshift(signals_b{pairs(i,1)}, tau_a(i));
    y =signals_b{pairs(i,2)};

    figure(1); subplot(size(pairs,1), 1, i);
    plot(x); hold on;  
    plot(y); 
    legend(sprintf('Signal %d', pairs(i,1)), sprintf('Signal %d', pairs(i,2))); 
    title ('Signals after the displacement')

    [acor, lag] = xcorr(x,y);
    lag_f = lag(:); %This is to avoid using the armonics
    acor_f = acor(850:1150);

    figure(2); subplot(size(pairs,1), 1, i);
    plot(lag_f,acor_f);
    title ('Autocorrelation of signals');

    [~, I] = max(abs(acor_f));
    tau_b(i) = -lag_f(I);

    travel_space(i) = 0.175*(2*(i-1)+1);
    travel_time(i) = (tau_b(i))/Fs_ref;
    velocity(i) = travel_space(i)/travel_time(i);
end
