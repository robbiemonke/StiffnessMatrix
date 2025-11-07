function [wrench] = getWrench(xi, params)
% This function takes an input parameters of the screw, xi (6x1), points,
% P (4x4), initial straight and cross lengths, ls_0 and lc_0, and the
% linear stiffness, k, and outputs the body wrench, wrench (6x1), of the
% system.
%
% Author: Robbie Monke
% Date: 7/28/2025

    % Parameters
    P = params.P;
    ls_0 = params.ls_0;
    lc_0 = params.lc_0;
    k = params.k;
    ls_max = params.ls_max;
    lc_max = params.lc_max;
    C1 = params.C1;
    C2 = params.C2;

    % Takes input of screw and converts to the body wrench
    T = screw2TMat(xi);

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
        
        % if lj < l_0
        %     fvec = zeros(3, 1);
        %     % disp("Slack")
        % else
        %     fvec = k * (lj - l_0) / (lj) * sj;
        %     % disp("Spring")
        % end

        fvec = k * (lj - l_0) / (lj ^ 2) * sj;

        F = F + fvec;
        M = M + vec2so3(rj) * (fvec);
    end
    wrench = [F; M];
end