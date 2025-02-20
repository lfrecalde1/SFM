function [H_norm] = normalization_values(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
x_mean = mean(X(1, :));
y_mean = mean(X(2, :));

x_var = var(X(1, :));
y_var = var(X(2, :));

% Compute normalization matrix
s_x = sqrt(2.0/x_var);
s_y = sqrt(2.0/y_var);

H_norm = [s_x, 0, -s_x*x_mean;...
          0, s_y, -s_y*y_mean;...
          0, 0, 1];
end