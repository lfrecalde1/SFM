function a_opt_vec = linefitCasadiL1norm(pts1, pts2, a_init)

addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;

% Convert to double precision if not already
A = double(pts1);
b = double(pts2);

% Decision variable: 9 parameters for the 3x3 H matrix
x = SX.sym('x', 2);    

%vector_difference = A*x - b;
cost = norm(A*x - b, 1);

% Define the NLP (no constraints g(x) in this basic example)
nlp = struct('x', x, 'f', cost);

% IPOPT solver options (you can tune these as needed)
opts = struct;
opts.print_time = 1;
opts.ipopt.print_level = 5;  % 0..12. Higher -> more verbose

% Create the CasADi solver using 'ipopt'
solver = nlpsol('solver', 'ipopt', nlp, opts);

x0 = a_init;
tic
% Solve
sol = solver('x0', x0);
toc
% Extract the solution
a_opt_vec = full(sol.x);  % Convert CasADi object to a numeric array
end