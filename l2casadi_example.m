clc, clear all, close all;

import casadi.*

% 1) Dimensions
n = 3;  % number of decision variables
m = 5;  % residual dimension (so r is in R^5)

% 2) Create an Opti instance

% 3) Decision variable
x = SX.sym('x', 3); 

% 4) Residual r = A*x - b (suppose we have A and b)
A = randn(m,n);
b = randn(m,1);
r = A*x - b;  % r will be an m-by-1 CasADi symbolic
cost = norm(r, 1);
% Define the NLP (no constraints g(x) in this basic example)
nlp = struct('x', x, 'f', cost);

% IPOPT solver options (you can tune these as needed)
opts = struct;
opts.print_time = 1;
opts.ipopt.print_level = 5;  % 0..12. Higher -> more verbose

% Create the CasADi solver using 'ipopt'
solver = nlpsol('solver', 'ipopt', nlp, opts);

x0 = [1;3;5];
tic
% Solve
sol = solver('x0', x0);
toc
% Extract the solution
a_opt_vec = full(sol.x);  % Convert CasADi object to a numeric array
a_opt_vec