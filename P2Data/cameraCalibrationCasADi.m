function [x_vector_opt, x_trans_opt, R_quaternion_opt, distortion_opt] = cameraCalibrationCasADi(pts1, pts2, A, x_init)

addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;

% Convert to double precision if not already
pts1 = double(pts1);
[size_optimization_x, size_optimization_y] = size(pts1);
pts2 = double(pts2);

% Decision variable: 9 parameters for the 3x3 H matrix
a_vector = SX.sym('full_estimation', 2 + (3+3) + (size_optimization_x-1)*size_optimization_y, 1);    

% Build the sum-of-squared-reprojection-errors cost
cost = 0;
%% Matrix of intrinsect parameters
%% Aux values for normalized
F_identity = [1, 0, 0;...
    0, 1, 0;...
    0, 0, 1];
Identity = [1, 0, 0, 0;...
    0, 1, 0, 0;...
    0, 0, 1, 0];

%% Reshape optimization variables
distortion = a_vector(1:2);
vector_optimization = a_vector(9:end);
x_vector = reshape(vector_optimization, 3, size_optimization_y);

%% Split the elements of the rotational and translational part
x_quaternion = a_vector(6:8);
x_trans = a_vector(3:5);

% 
% 
%% Section to map the R^n vector to quaternions
quaternion = [(1 - x_quaternion(1:3)'*x_quaternion(1:3))/(1  + x_quaternion(1:3)'*x_quaternion(1:3)); 2*x_quaternion(1)/(1 + x_quaternion(1:3)'*x_quaternion(1:3)); 2*x_quaternion(2)/(1 + x_quaternion(1:3)'*x_quaternion(1:3)); 2*x_quaternion(3)/(1 + x_quaternion(1:3)'*x_quaternion(1:3))];
trans_aux = [x_trans(1); x_trans(2); x_trans(3)];
% 
% REAL VALUES
U_real = pts2(1:2, :);

% TRANSFORMATION
T_estimated = [quatTorot(quaternion), trans_aux; 0 0 0 1];
% 
% NORMALIZED VALUES
values_normalized = F_identity*Identity*T_estimated*[x_vector; ones(1, size(pts1, 2))];
aux_normalization = [values_normalized(3, :); values_normalized(3, :)];
values_normalized_aux = values_normalized(1:2, :)./aux_normalization;

%     % APLYING BARREL DISTORTION
radius = sqrt(sum(values_normalized_aux.^2, 1));
D = 1 + distortion(1)*radius.^2 + distortion(2)*radius.^4;
D_aux = [D;D];
x_warp = values_normalized_aux.*D_aux;
x_warp_aux = [x_warp; ones(1, length(x_warp))];

% 
U_improved = A * x_warp_aux;
U_normalized_aux = [U_improved(3,:); U_improved(3,:)];
U_improved_final = U_improved(1:2,:) ./ U_normalized_aux;

error = (U_real - U_improved_final);

error_reshape = reshape(error, 2*size(error,2), 1);

cost = cost + error_reshape'*error_reshape;

% Define the NLP (no constraints g(x) in this basic example)
nlp = struct('x', a_vector, 'f', cost);

% IPOPT solver options (you can tune these as needed)
opts = struct;
opts.print_time = 1;
opts.ipopt.print_level = 5;  % 0..12. Higher -> more verbose

% Create the CasADi solver using 'ipopt'
solver = nlpsol('solver', 'ipopt', nlp, opts);
% 
x0 = x_init;
% % 
% Solve
sol = solver('x0', x0);
% 
% Extract the solution
a_opt_vec = full(sol.x);  % Convert CasADi object to a numeric array
% Getting vlues from solution
distortion_opt = a_opt_vec(1:2);
vector_optimization_opt = a_opt_vec(9:end);
x_vector_opt = reshape(vector_optimization_opt, 3, size_optimization_y);

%% Split the elements of the rotational and translational part
x_quaternion_opt = a_opt_vec(6:8);
quaternion_opt = [(1 - x_quaternion_opt(1:3)'*x_quaternion_opt(1:3))/(1  + x_quaternion_opt(1:3)'*x_quaternion_opt(1:3)); 2*x_quaternion_opt(1)/(1 + x_quaternion_opt(1:3)'*x_quaternion_opt(1:3)); 2*x_quaternion_opt(2)/(1 + x_quaternion_opt(1:3)'*x_quaternion_opt(1:3)); 2*x_quaternion_opt(3)/(1 + x_quaternion_opt(1:3)'*x_quaternion_opt(1:3))];

x_trans_opt = a_opt_vec(3:5);
R_quaternion_opt = quatTorot(quaternion_opt);

end