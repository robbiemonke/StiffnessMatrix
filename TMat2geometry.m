function [u, theta, rho, h] = TMat2geometry(T)
% This function takes argument of a transformation matrix (4x4) and outputs
% the geometrical parameters of the axis of rotation u (3x1), angle of
% rotation theta, instantaneous center of rotation rho, and pitch of the
% screw h.
%
% Author: Robbie Monke
% Date: 4/13/2025
R = T(1:3, 1:3);
t = T(1:3, 4);

if norm(R - eye(3)) < 1e-12
    theta = 0;
    u = t / norm(t);
    h = norm(t);
    rho = NaN;
else
    [u, theta] = EquivalentAxis3(R);
    rho = (eye(3) - R') * t / (2 * (1 - cos(theta)));
    h = u' * t / theta;
end

end