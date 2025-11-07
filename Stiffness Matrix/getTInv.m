function [Tinv] = getTInv(T)
% This function takes an input of a rotation matrix, R (3x3), and outputs
% the equivalent axis of rotation, u (3x1), and angle of rotation, theta.
% 
% Author: Robbie Monke
% Date: 2/9/2025
R = T(1:3, 1:3);
t = T(1:3, 4);
Tinv = [R' -R'*t; 0 0 0 1];
end