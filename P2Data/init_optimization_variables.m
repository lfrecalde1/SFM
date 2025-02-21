function [X_init] = init_optimization_variables(translation, rotation, points)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X_init = [];
quaternion_estimated = rotm2quat(rotation);
x_quaternion = quaternion_estimated(2:4)/quaternion_estimated(1);
X_init = [X_init, 0.0, 0.0, translation', x_quaternion];
for k=1:size(points, 2)
    X_init = [X_init,points(1:3, k)'];
end
end