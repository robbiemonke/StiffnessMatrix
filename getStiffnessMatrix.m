function [K] = getStiffnessMatrix(xi, P, ls_0, lc_0, k)
% This function takes an input parameters of the screw, xi (6x1), points,
% P (4x4), initial straight and cross lengths, ls_0 and lc_0, and the
% linear stiffness, k, and outputs the stiffness matrix, K (6x6), of the
% system.
%
% Author: Robbie Monke
% Date: 7/28/2025
    T = screw2TMat(xi);

    % Connection Matrices for original design A1A2 D1D2
    C1 = [1 1 1 1 0 0 0 0 0 0 0 0;...
          0 0 0 0 1 1 0 0 0 0 0 0;...
          0 0 0 0 0 0 1 1 0 0 0 0;...
          0 0 0 0 0 0 0 0 1 1 1 1];
    
    C2 = [1 0 0 0 1 0 1 0 1 0 0 0;...
          0 1 0 0 0 0 0 0 0 1 0 0;...
          0 0 1 0 0 0 0 0 0 0 1 0;...
          0 0 0 1 0 1 0 1 0 0 0 1];
    
    
    P1 = P * C1;
    P2 = P * C2;
    S = P1 - T * P2;
    
    dT_dxi = screw2dT(xi);
    for j = 1:12
        sj = S(1:3, j);
        lj = norm(sj);
        pj = P2(:, j);
        
        % Determination of initial length depending on if its a cross member or
        % straight member.
        if j == 1 || j == 4 || j == 9 || j == 12
            l_0 = lc_0;
        else
            l_0 = ls_0;
        end

        % Partial force partial s
        dfvec_ds = k * (((lj - l_0) / lj ^ 2) * eye(3) + ((2 * l_0 - lj) / lj ^ 4) * (sj * sj'));
    
        % Partial derivative of s w.r.t. xi
        ds_dxi = zeros(3, 6);
        for i = 1:6
            ds_dxi(:, i) = -dT_dxi(1:3, :, i) * pj;
        end
    
        % Development of moment arm skew symmetric matrix
        rj = P1(1:3, j);
        rj_hat = vec2so3(rj);
        A = [eye(3); rj_hat];
    
        Kj(:, :, j) = A * dfvec_ds * ds_dxi;
    end
    
    K = sum(Kj, 3);
end