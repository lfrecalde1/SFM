function [pixels_aux_estimated] = projection_values(H1, pts3D_4xN, k1, k2, K)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F_identity = [1, 0, 0;...
              0, 1, 0;...
              0, 0, 1];
Identity = [1, 0, 0, 0;...
            0, 1, 0, 0;...
            0, 0, 1, 0];

values_normalized_estimated = F_identity*Identity*H1*pts3D_4xN;
values_normalized_estimated = values_normalized_estimated(1:2, :)./values_normalized_estimated(3, :);
radius_estimated = vecnorm(values_normalized_estimated);
D_estimated = 1 + k1*radius_estimated.^2 + k2*radius_estimated.^4;
x_warp_estimated = values_normalized_estimated.*D_estimated;
x_warp_estimated = [x_warp_estimated; ones(1, length(x_warp_estimated))];
pixels_aux_estimated = K * x_warp_estimated;
pixels_aux_estimated = pixels_aux_estimated(1:2,:) ./ pixels_aux_estimated(3,:);

end