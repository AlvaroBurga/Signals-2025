%% PDF de X = cos(theta) con theta ~ Unif(0, 2pi)
clear; clc; close all;

% Número de muestras
N = 1e6;

% Generar theta uniforme
theta = 2*pi*rand(N,1);

% Variable aleatoria X = cos(theta)
X = cos(theta);

%% Estadísticos (deben dar mu=0, var=1/2)
muX = mean(X);
varX = var(X);

%% Rango para PDFs
x_vals = linspace(-1, 1, 1000);

%% PDF teórica de X (arcsine distribution)
f_arcsin = 1 ./ (pi * sqrt(1 - x_vals.^2));
% Nota: esto explota en los bordes, Matlab lo maneja

%% Normal con misma media y varianza
sigmaX = sqrt(varX);
f_normal = (1./(sigmaX*sqrt(2*pi))) .* exp(-(x_vals - muX).^2 ./ (2*sigmaX^2));

%% Graficar histogramas + pdfs
figure;
hold on;

% Histograma normalizado
histogram(X, 200, 'Normalization', 'pdf', 'DisplayStyle', 'bar', ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');

% PDF teórica arco-seno
plot(x_vals, f_arcsin, 'r', 'LineWidth', 2);

% PDF normal equivalente
plot(x_vals, f_normal, 'k--', 'LineWidth', 2);

xlabel('x');
ylabel('PDF');
title('PDF de X = cos(\theta) y Normal con misma media/varianza');
legend('Histograma de X', 'PDF Arcsine teórica', 'Normal equivalente');
grid on;
