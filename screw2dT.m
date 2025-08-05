function [dT_dxi] = screw2dT(xi)
% This function takes an input of the body screw, xi (6x1) and calculates the
% derivative of the transformation matrix w.r.t. the screw. The output
% is stacked 4x4 matrices of each partial derivative resulting in a 4x4x6 
% output, dT_dxi.
%
% Author: Robbie Monke
% Date: 7/24/2025
    % Split of given screw into angular and linear components
    omega = xi(1:3);
    v = xi(4:6);
    omega_hat = vec2so3(omega);
    theta = norm(omega);
    
    % Rotation matrix and its derivative w.r.t. its exponential coordinates
    dR = omega2dR(omega);

    % Definition of coefficients for left jacobian of SO(3)
    % Luo et al. - 2022 - The Geometry and Kinematics of the Matrix Lie Group SE(3)
    % Page 8 Eq (15)

    if theta < eps % Pure translation case
        a = eye(3);
        b = 1 / 2;
        c = 1 / 6;
        da_dtheta = zeros(3, 3);
        db_dtheta = 0;
        dc_dtheta = 0;
    else % Translation and rotation case
        a = eye(3);
        b = ((1 - cos(theta)) / theta ^ 2);
        c = ((theta - sin(theta)) / theta ^ 3);
        da_dtheta = zeros(3, 3);
        db_dtheta = (theta * sin(theta) + 2 * cos(theta) - 2) / theta ^ 3;
        dc_dtheta = (3 * sin(theta) - theta * cos(theta) - 2 * theta) / theta ^ 4;
    end

    G = a + b * omega_hat + c * omega_hat ^ 2;

    for j = 1:3
        % Derivative of omega_hat
        e_j = zeros(3,1); 
        e_j(j) = 1;
        domegahat_domega = vec2so3(e_j);

        if theta < eps % Pure translation case
            % Derivative of theta w.r.t. omega j
            dtheta_domegaj = 0;

            % Calculation of derivative of the left jacobian w.r.t. omega j
            dJ_domega = b * domegahat_domega + ...
                c * (domegahat_domega * omega_hat + omega_hat * domegahat_domega);

            % Forcing values of zeros for pure translation case
            dt_domega = zeros(3, 1);
            
        else % Translation and rotation case
            % Derivative of theta w.r.t. omega j
            dtheta_domegaj = omega(j) / (theta);

            % Calculation of derivative of the left jacobian w.r.t. omega j
            dJ_domega = da_dtheta * dtheta_domegaj * eye(3) + ...
                db_dtheta * dtheta_domegaj * omega_hat + ...
                b * domegahat_domega + ...
                dc_dtheta * dtheta_domegaj * omega_hat * omega_hat + ...
                c * (domegahat_domega * omega_hat + omega_hat * domegahat_domega);

            % Multiplication of derivative of left jacobian and v to calculate
            % the derivative of the translation t w.r.t. omega
            dt_domega = dJ_domega * v;
        end

        % Partial of T w.r.t. omega (4x4x3)
        dT_domega(:, :, j) = [dR(:, :, j), dt_domega; zeros(1, 4)];

        % Partial T partial v (4x4x3)
        dT_dv(:, :, j) = [zeros(3, 3) G(:, j); zeros(1, 4)];
    end

    % Concatenation of stacked matrices (4x4x6)
    dT_dxi = cat(3, dT_domega, dT_dv);
end