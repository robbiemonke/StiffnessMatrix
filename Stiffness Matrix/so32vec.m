function v = so32vec(S)
% Converts a skew-symmetric matrix S (3x3) into a vector v (3x1)
% S is an input of a 3x3 skew symmetric matrix.
% v is an output of a 3x1 vector.
% 
% Author: Robbie Monke
% Date: 2/9/2025

v = [S(3, 2); S(1, 3); S(2, 1)];
end
