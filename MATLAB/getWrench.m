function [wrench] = getWrench(xi, P, ls_0, lc_0, k)
% This function takes an input parameters of the screw, xi (6x1), points,
% P (4x4), initial straight and cross lengths, ls_0 and lc_0, and the
% linear stiffness, k, and outputs the body wrench, wrench (6x1), of the
% system.
%
% Author: Robbie Monke
% Date: 7/28/2025
    % Takes input of screw and converts to the body wrench
    T = screw2TMat(xi);

    % Definition of connection matrices C1, C2 for A1A2 D1D2 configuration
    C1 = [1 1 1 1 0 0 0 0 0 0 0 0;...
          0 0 0 0 1 1 0 0 0 0 0 0;...
          0 0 0 0 0 0 1 1 0 0 0 0;...
          0 0 0 0 0 0 0 0 1 1 1 1];

    C2 = [1 0 0 0 1 0 1 0 1 0 0 0;...
          0 1 0 0 0 0 0 0 0 1 0 0;...
          0 0 1 0 0 0 0 0 0 0 1 0;...
          0 0 0 1 0 1 0 1 0 0 0 1];

    % Solves equation for connection members vectors
    S = P * C1 - T * P * C2;

    F = [0 0 0]';
    M = [0 0 0]';
    for j = 1:12
        sj = S(1:3, j);
        lj = norm(sj);
        R = P * C1;
        rj = R(1:3, j);

        if j == 1 || j == 4 || j == 9 || j == 12
            l_0 = lc_0;
        else
            l_0 = ls_0;
        end
        
        fvec = k * (lj - l_0) / (lj ^ 2) * sj;
        F = F + fvec;
        M = M + vec2so3(rj) * (fvec);
    end
    wrench = [F; M];
end
