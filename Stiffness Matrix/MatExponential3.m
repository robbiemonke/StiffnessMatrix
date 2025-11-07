function [R] = MatExponential3(omega, theta)
% This function takes parameters of a 3x1 matrix omega representing the
% axis of rotation and an angle of rotation, theta, and outputs the equivalent
% rotation matrix R.
% 
% The angle theta is given in units of radians.
% 
% Author: Robbie Monke
% Date: 2/6/2025

% Calculation of the rotation matrix (R) from the axis of rotation (omega) and
% angle of rotation (theta).
if norm(omega) == 0
    R = eye(3);
else
    omega = omega / norm(omega);
    R = eye(3) + sin(theta) * vec2so3(omega) + (1 - cos(theta)) * (vec2so3(omega) ^ 2);
end
end