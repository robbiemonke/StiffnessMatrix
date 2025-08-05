function [R] = MatExponential3(omega, theta)
% This function takes parameters of a 3x1 matrix omega representing the
% axis of rotation and an angle of rotation, theta, and outputs the equivalent
% rotation matrix R.
% 
% The angle theta is given in units of radians.
% 
% Author: Robbie Monke
% Date: 2/6/2025

% Conversion of axis of rotation to a unit vector
omega = omega/norm(omega);

% Creation of omega_hat; The so3 representation of the vector omega
omega_hat = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];

% Calculation of the rotation matrix (R) from the axis of rotation (omega) and
% angle of rotation (theta).
R = eye(3) + sin(theta) * omega_hat + (1 - cos(theta)) * omega_hat * omega_hat;
end