clc, clear all, close all;

import casadi.*

% 1) Dimensions
n = 3;  % number of decision variables
m = 5;  % residual dimension (so r is in R^5)

% 2) Create an Opti instance
opti = Opti();

% 3) Decision variable
x = opti.variable(n,1);

% 4) Residual r = A*x - b (suppose we have A and b)
A = randn(m,n);
b = randn(m,1);
r = A*x - b;  % r will be an m-by-1 CasADi symbolic

% 5) Slack variables for L1 absolute value
s = opti.variable(m,1);  % s_i >= 0

% 6) Impose constraints to enforce s_i >= |r_i|
opti.subject_to( r <= s );
opti.subject_to(-r <= s );
opti.subject_to( s >= 0 );

% 7) Objective: sum of slack variables
obj = sum(s);
opti.minimize(obj);

% Additional linear constraints? Add them similarly:
% e.g., x >= 0 => x(1) >= 0, x(2) >=0, etc:
% opti.subject_to(x >= 0);

% 8) Choose a solver
% (For an LP, you might prefer a dedicated LP/QP solver; 
%  Ipopt can still handle it as an NLP, but it's not specialized for pure LP.)
opti.solver('ipopt');

% 9) Solve
tic
sol = opti.solve();
toc
x_opt = sol.value(x);
s_opt = sol.value(s);

disp('Optimal x:');
disp(x_opt);
disp('Optimal cost (sum of absolute residuals):');
disp(sum(s_opt));