%% SLP 2025/2026 - Excercize 1

clc; clear;close all;

%% Ex1

fo = 1e3; %Hz - Carrier frequency
Fs = 20*fo; %Hz - Sampling Frecuency
t = 0:1/Fs:0.1;

xn = sin(2*pi*fo*t);

figure;
plot(t, xn, 'LineWidth',2); %ALWAYS USE THE LABELS!!
xlabel('Time (s)');
ylabel('Amplitude');
title('Sine Wave at 1 kHz');
grid on;

%% Vibration
clc; clear;close all;

load("Experiment_table.mat")

plot(x3(:,1), x3(:,4)); hold on
plot(x4(:,1), x4(:,4));
xlabel('Time [s]');
ylabel('Amplitude');

xn = x1(:, 2); %This is the emmited signal
t = x1(:,1); %The time of the emmited signal
Ts = t(3) - t(2); %Sampling period
Fs = 1/Ts; %Sampling frecuency
N = length(xn); % Number of samples


Xf = fftshift(fft(xn))/sqrt(N); %% Muy importante poner el fftship para tener el espectro centrado. Igualmente normalizar dividiento por sqrt(N)
fn = (-N/2 + 1 : N/2)/N; % Normalized frequency vector (va de -0.5 a 0.5) 
f = fn * Fs; % Real frecuency vector

figure
plot(f, abs(Xf), 'LineWidth',2)
xlabel('Frecuency (Hz)');
ylabel('Amplitude');

% define frecuency response
Hf = ones(N, 1); % Filter to remove the DC component
Hf(f > -5 & f<5) = 0; % Our signal is in 11Hhz, so we filter what we have beforeç

figure; hold on
plot(f, abs(Hf), 'LineWidth',2);
plot(f, abs(Xf), 'LineWidth',2)
xlabel('Normalize Frecuency (Hz)');
ylabel('Amplitude');

% ifft:
hn = ifft(ifftshift(Hf)); % Nunca olvidar el shift (Señal del filtro en el tiempo)

figure;
plot(t, hn, 'LineWidth',2)
xlabel('Time (s)');
ylabel('Amplitude');

hn_t = hn; %Optimized response reduce the computational workload
hn_t(abs(hn) < 0.01) = 0;
hold on;
plot(t, hn_t)

%Compute convolution (apply the filter)
yn = conv(xn, hn); %Filtered response (without the DC component)
L = N + N -1; %Length of the convolution matrix
tc = (0: L-1) * Ts; %Time vector for the convolution matrix

% Convolution by con matrix
Hg = zeros(L, N); %Convolution matrix of the system
Xg = zeros(L, N); %Convolution matrix of the emmited signal
for i = 1 : N
    Hg((1:N)+ i-1 , i) = hn;
    Xg((1:N)+ i-1 , i) = xn;

end

ynm = Hg * xn; %Getting the response via convolution matrix

figure;
plot(tc, yn, 'LineWidth',2)
plot(tc, ynm, 'LineWidth',2, "Color","r")
xlabel('Time (s)');
ylabel('Amplitude');

% deconvolution
%Here we are creating the inverse matrix of convolution
Hinverse = (Hg' * Hg) \ Hg'; % Use '\' apply to ^(-1) of a matix. The ' means the Hermitian (Transpose of a complex number) = Tranpose + conj.
xn_est =  Hinverse * yn; 
figure; hold on;
plot(t, xn, '-o', 'LineWidth',2)
plot(t, xn_est, '-x', 'LineWidth',2, "Color","r")
xlabel('Time (s)');
ylabel('x[n]');


% System iDENTIFICATION
%Here we are creating the inverse matrix of convolution
Xinverse = (Xg' * Xg) \ Xg'; % Use '\' apply to ^(-1) of a matix. The ' means the Hermitian (Transpose of a complex number) = Tranpose + conj.
hn_est =  Xinverse * yn; 
figure; hold on;
plot(t, hn, '-o', 'LineWidth',2)
plot(t, hn_est, '-x', 'LineWidth',2, "Color","r")
xlabel('Time (s)');
ylabel('x[n]');

%% Optimiation and Quadratic Form
clc; clear;close all;

x = [3, 4]';
h = [1, .5, .2, -2, 3, -.1];

N = length(x);
M = length(h);
L = N + M -1;
H = zeros(L, N);
for i = 1 : N
    H((1:M) + i - 1, i) = h;
end


% Least square Fitting
x1 = -10 : .5 : 10;
x2 = -10 : .5 : 10;
K1 = length(x1);
K2 = length(x2);

%Alternative 1
w = randn(L, 1)*0.1; %Your friend, the noise :) (Or maybe the enemy)
y = H * x + w;
mse = zeros(K1, K2);

for k1 = 1 : K1
    for k2 = 1 : K2
        % Perform least squares fitting for each combination of x1 and x2
        x_trial = [x1(k1), x2(k2)]';
        mse(k1, k2) = mean(abs(y - H * x_trial).^2);
    end
end

figure;
surf(x1, x2, mse)

%Montecarlo simulation
for mc = 1:1e5
    w = randn(L, 1)*0.1; %Your friend, the noise :) (Or maybe the enemy)
    y = H * x + w;
    x_est(mc, :) = (H' * H) \ H' * y;
end

figure; hold on;
plot(x_est(:,1), x_est(:,2), '.')
plot(x(1), x(2), 'r*', 'MarkerSize', 5)
xlabel('x1')
ylabel('x2')
axis equal

meanXest = mean (x_est());
error = meanXest-x'

%% Experiment
xn1 = x1(:,2);
t1 = x1(:,1);

xn2 = x2(:,2);
t2 = x2(:,1);

figure; hold on;
plot(t1, xn1);
plot(t2, xn2);
xlabel('time [s]')

N1 = length(x1);
N2 = length(x2);
L = N1 + N2 -1;
tau = (-N1 + 1 : N1 -1) * Ts;
Rx1x2 = xcorr(xn1, xn2); %With corelation you solve non-linear equations

figure;
plot(tau, abs(Rx1x2));
xlabel('Delay [s]')