% Calculates the Transformation Matrix given a twist - screw (omega, v) and
% its magnitude (theta)
%
% omega: 3x1 vector (rotational component of the screw) - unit vector
% theta: 1x1 scalar magrinude of the screw
% v: 3x1 vector (translational velocity component of the screw)
% 
% T: 4x4 transformation matrix
% 
% 
function [T] = TMatExponential3(u,theta,vtilde) 
    zeroThreshold = 1e-8;
    if and(norm(u) < zeroThreshold,isnumeric(u))
        R = eye(3);
        p = vtilde * theta;
    else
        R = MatExponential3(u,theta);
        p = (eye(3) - R) * vec2so3(u) * vtilde + u * u' * vtilde * theta;
    end
    T = [R, p; zeros(1,3), 1];
end
        