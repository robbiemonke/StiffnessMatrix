function [T] = geometry2TMat(u, theta, rho, h)
% This function takes a parameter of a 1x3 direction of translation (u),
% an angle of rotation (theta), a 1x3 axis of rotation (rho), and a pitch (h),
% and outputs the 4x4 translation matrix (T).
% 
% Author: Robbie Monke
% Date: 3/15/2025
R = MatExponential3(u, theta);
t = (eye(3)-R)*rho' + h*theta*u';
T = [R, t; 0 0 0 1];
end
