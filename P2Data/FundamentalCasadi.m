function [F_opt] = FundamentalCasadi(X1, X2, x_init)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');

import casadi.*


x = SX.sym('x', 9);
F = [];
for k=1:size(X1, 2)

    aux_f = kron(X1(1:3, k), X2(1:3, k))';
    F = [F; aux_f];
end
r = F*x;
cost = r'*r;
g = x'*x;

% Define the NLP (no constraints g(x) in this basic example)
nlp = struct('x', x, 'f', cost, 'g', g);

% IPOPT solver options (you can tune these as needed)
opts = struct;
opts.print_time = 0;
opts.ipopt.print_level = 5;  % 0..12. Higher -> more verbose

% Create the CasADi solver using 'ipopt'
solver = nlpsol('solver', 'ipopt', nlp, opts);
x0 = reshape(x_init, 9, 1);

tic
% Solve
sol = solver('x0', x0, 'lbg', 1, 'ubg',1);

toc
F_opt_vec = full(sol.x);  % Convert CasADi object to a numeric array
F_opt = reshape(F_opt_vec, 3, 3);
end