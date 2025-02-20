function [F_final] = fundamental_analytical(X, U, K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Normalize the data
[H_X] = normalization_values_new();
[H_U] = normalization_values_new();
X = inv(K) * X;
U = inv(K) * U;

F = [];
for k=1:size(U, 2)

    aux_f = kron(X(1:3, k), U(1:3, k))';
    F = [F; aux_f];
end

% Perform Singular Value Decomposition
[U_f, S_f, V_f] = svd(F);

% Extract singular values from the diagonal of S into a vector
sing_vals = diag(S_f);

% Find the index of the smallest singular value
[~, idx] = min(sing_vals);

% The solution is the corresponding column of V
F_nomr = V_f(:, idx);
F_norm = reshape(F_nomr, 3, 3)';
F_final = pinv(K)*F_norm*K;
%F_final = F_norm;
end