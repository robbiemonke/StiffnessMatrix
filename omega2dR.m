function [dR] = omega2dR(omega)
% This function calculates the derivatives of a rotation, R(omega) =
% exp(vec2so3(omega)), with respect to its exponential coordinates omega.
% 
% This function stems from equation (7) of http://dx.doi.org/10.1007/s10851-014-0528-x
%
% omega - angular velocity 3x1
% dR - derivatives of the rotation 3x3x3
%
% Author: Robbie Monke
% Date: 7/23/2025
    theta = norm(omega);
    omega_bar = omega / theta; % Normalization of omega for calculation

    if theta < eps % Pure translation case
        dR = zeros(3, 3, 3);
        return;
    end

    % Calculation of 3 layers of derivative of the rotation matrix w.r.t.
    % its exponential coordinates
    dR = zeros(3, 3, 3); % Pre-allocation
    for i = 1:3
        ei = zeros(3, 1);
        ei(i) = 1;

        dR(:, :, i) = cos(theta) * omega_bar(i) * vec2so3(omega_bar) + ...
             sin(theta) * omega_bar(i) * (vec2so3(omega_bar) ^ 2) + ...
             (sin(theta) / theta) * vec2so3(ei - omega_bar(i) * omega_bar) + ...
             ((1 - cos(theta)) / theta) * (ei * omega_bar' + omega_bar * ei' - 2 * omega_bar(i) * (omega_bar * omega_bar'));
    end
end