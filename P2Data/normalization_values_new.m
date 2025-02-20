function [H_norm] = normalization_values_new()
%NORMALIZATION_VALUES  Normalizes 2D image coordinates to [-1,1].

% The image center
cx = 640;   % Half of 1280
cy = 512;   % Half of 1024

% Scale so that x ranges from -1 to +1 and y ranges from -1 to +1
% x' = (x - cx) / cx  ==> x' = (1/cx)*x - (cx/cx) = (1/cx)*x - 1
% y' = (y - cy) / cy  ==> y' = (1/cy)*y - 1
sx = 1 / cx;   % 1/640
sy = 1 / cy;   % 1/512

% Build the normalization matrix in homogeneous form
H_norm = [ sx   0   -1;    % x' = sx*x + 0*y + (-1)
           0    sy  -1;    % y' = 0*x + sy*y + (-1)
           0    0    1 ];

end