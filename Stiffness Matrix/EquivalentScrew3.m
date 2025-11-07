% This function takes a parameters of the 4x4 transformation matrix (T) and
% outputs the 3x1 axis of rotation (u), 3x1 vtilde, and angle of rotation
% (theta).
% Date: 3/15/2025
function [omega,theta,v] = EquivalentScrew3(T)
    R = T(1:3,1:3); p = T(1:3,4);
    if trace(R)==3 % No rotation
        omega =[0,0,0]';
        theta = norm(p);
        v = p/theta;
    else
    [omega,theta] = EquivalentAxis3(R);
    omegabracket = vec2so3(omega);
    Ginv =(eye(3)/theta - (omegabracket/2)+ (1/theta-(cot(theta/2)/2))*omegabracket^2);
    v = Ginv*p;
    end
end