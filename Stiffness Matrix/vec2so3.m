function S = vec2so3(v)
% Converts a vector v (3x1) into a skew-symmetric matrix S (3x3)
% v is an input of a 3x1 vector.
% S is an output of a 3x3 skew symmetric matrix.
% 
% Author: Robbie Monke
% Date: 2/9/2025

S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
