function [E_total, S, C1, C2] = fFindFunction(xi, params, Fext)
    
    % Parameters
    P = params.P;
    ls_0 = params.ls_0;
    lc_0 = params.lc_0;
    k = params.k;
    C1 = params.C1;
    C2 = params.C2;

    % Takes input of screw and converts to the transformation matrix
    T = screw2TMat(xi);

    % Solves equation for connection members vectors
    S = P * C1 - T * P * C2;
    S3 = S(1:3, :);

    % Iterates through connection members and solves for the total energy
    % in the system
    E_total = 0;
    E = zeros(length(S), 1);
    for i = 1:length(S)
        if i == 1 || i == 4 || i == 9 || i == 12 % Member is diagonal
            if (norm(S3(:, i)) - lc_0) <= 0 % Member is slack
                E(i) = 0;
            else % Member is in tension
                E(i) = 0.5 * k * (norm(S3(:,i)) - lc_0) ^ 2;
            end
        else % Member is straight
            if (norm(S3(:, i)) - ls_0) <= 0 % Member is slack
                E(i) = 0;
            else % Member is in tension
                E(i) = 0.5 * k * (norm(S3(:,i)) - ls_0) ^ 2;
            end
        end
        E_total = E(i) + E_total;
    end

    F = [Fext; vec2so3(P(1:3, 1)) * Fext];
    % F = [0.01 0 0 0 0 0]';
    E_total = E_total + xi' * F;
end