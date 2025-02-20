%% Code to theck the behavior of ransac to outliers
clc, clear all, close all;


% True parameters for the linear model
a_true = 2.5;
b_true = 1.0;

% Generate x-values
num_points = 30;
x = linspace(0, 10, num_points);


% Calculate y-values
y = a_true .* x + b_true;


%% Introduce outliers
rng('default');  % For reproducible results

% Increase this to get more outliers
num_outliers = 5;  
outlier_indices = randperm(num_points, num_outliers);

% Increase or decrease these points by a large random amount
% Adjust the multiplier (e.g., 25, 50, etc.) to control severity of outliers
multiplier = 200; 
y(outlier_indices) = y(outlier_indices) + multiplier * randn(size(outlier_indices));


%% Create data for matrix for analytical optimization
A = [x', ones(size(x,2), 1)];
Y = [y'];
parameters = pinv(A)*Y

%% Computing optimal values that fit the line 
x_opti = linefitCasadiL1norm(A, Y, [1; 1])

tic
%% Computing parameters based on ransac
best_model = fit(A, Y, 2, std(Y)/2)
toc
% Plot the data
figure;
plot(x, y, 'bo', 'MarkerFaceColor','b');  % Plot data points
hold on;

% Plot the true line (without noise)
x_fit = linspace(0, 10, 100);    % more points for a smooth line
y_fit = parameters(1) .* x_fit + parameters(2);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Plot the true line (without noise)
x_fit_casadi = linspace(0, 10, 100);    % more points for a smooth line
y_fit_casadi = x_opti(1) .* x_fit + x_opti(2);
plot(x_fit_casadi, y_fit_casadi, 'g-', 'LineWidth', 2);

% Plot the true line (without noise)
x_fit_ransac = linspace(0, 10, 100);    % more points for a smooth line
y_fit_ransac = best_model(1) .* x_fit + best_model(2);
plot(x_fit_ransac, y_fit_ransac, 'y--', 'LineWidth', 2);

xlabel('x');
ylabel('y');
legend('Noisy Data', 'True Line');
title('Synthetic Linear Data with Noise');
grid on;
hold off;