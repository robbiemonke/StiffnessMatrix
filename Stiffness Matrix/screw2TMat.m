function [T] = screw2TMat(S)
    epsilon = 1e-8;
    if norm(S(1:3))<epsilon % omega =0
        omega   = zeros(3,1);
        theta   = norm(S(4:6));
        v       = S(4:6)/theta;
    else
        theta   = norm(S(1:3));
        omega   = S(1:3)/theta;
        v       = S(4:6)/theta;
    end
    T = TMatExponential3(omega,theta,v);
end