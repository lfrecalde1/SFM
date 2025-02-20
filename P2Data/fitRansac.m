function [best_model, iterations_done] = fitRansac(X1, X2, num_sample, F_init, threshold)

% Initialization
num_iterations  = Inf;
iterations_done = 0;

max_inlier_count = 0;
best_model       = [];

desired_prob  = 0.95;

F = [];
Y = [];
for k=1:size(X1, 2)
    aux_f = kron(X1(1:3, k), X2(1:3, k))';
    F = [F; aux_f];

    aux_y = 0;
    Y = [Y; aux_y];
end
% threshold = std(Y)/2;
% Combine A and Y into one data matrix for easier row shuffling
total_data = [F, Y];   % size: (N x (d+1))
data_size  = size(total_data, 1);


% RANSAC loop (adaptive number of iterations)
while num_iterations > iterations_done
    % Shuffle the rows
    idx = randperm(data_size);
    total_data = total_data(idx, :);
    

    % Take the first 'num_sample' as the subset
    sample_data = total_data(1:num_sample, :);

    % Fit model to the sample subset
    F_samp = sample_data(:, 1:end-1);
    Y_samp = sample_data(:, end);
    
    [F_optimization] = FundamentalCasadiAux(F_samp, Y_samp, F_init);
    % 
    % % Compute error across all data
    % % (Ensure dimensions align so that A * model is (N x 1))
    estimated_model = reshape(F_optimization, 9, 1);
    y_estimated = F * estimated_model;   % predicted Y
    err   = abs(Y - y_estimated);
    % 
    % % Count how many points are inliers (within threshold)
    inlier_count = nnz(err < threshold);

    % Update best model if this one has more inliers
    if inlier_count > max_inlier_count
        max_inlier_count = inlier_count;
        best_model       = F_optimization;
    end
    % 
    % % Recompute probability of outlier and number of iterations needed
    prob_outlier  = 1 - (inlier_count / data_size);
    num_iterations = log(1 - desired_prob) / ...
        log(1 - (1 - prob_outlier)^num_sample);
    
    iterations_done = iterations_done + 1;
    F_init = F_optimization;
end

end