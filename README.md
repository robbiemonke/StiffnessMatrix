# StiffnessMatrix
Repository of current work for stiffness matrix derivation of tensegrity systems. Includes the MATLAB code using and testing the derivation.

# MATLAB Code:
omega2dR.m - derivative of rotation matrix w.r.t. its exponential coordinates.

screw2dT.m - derivative of transformation matrix w.r.t. its exponential coordinates.

TestScrew2dT.mlx - Test script comparing screw2dT.m against finite difference method.

getStiffnessMatrix.m - calculation of stiffness matrix using derived method.

TestStiffnessMatrix.mlx - Test script comparing getStiffnessMatrix.m against finite difference method.
