% Constants and Params
% All units in kg, cm, and s
% Connection Matrices
C1 = [1 1 1 1 0 0 0 0 0 0 0 0;...
      0 0 0 0 1 1 0 0 0 0 0 0;...
      0 0 0 0 0 0 1 1 0 0 0 0;...
      0 0 0 0 0 0 0 0 1 1 1 1];

C2 = [1 0 0 0 1 0 1 0 1 0 0 0;...
      0 1 0 0 0 0 0 0 0 1 0 0;...
      0 0 1 0 0 0 0 0 0 0 1 0;...
      0 0 0 1 0 1 0 1 0 0 0 1];

% First arc's points (cm)
% wrt origin
P = [-3.6195, 0, 0, 1;...
    -2.2924, -2.4130, 0, 1;...
    2.2924, -2.4130, 0, 1;...
    3.6195, 0, 0, 1]';

% Straight tendon member length (cm)
ls_0 = 3.175;
ls_max = ls_0 + 2;

% Diagonal tendon member length (cm)
lc_0 = 3.810;
lc_max = lc_0 + 2;

% Spring constant (kg cm / cm s^2)
k = 0.77 * 100;

% Spatial Intertia Matrix (kg cm ^ 2)
m = 23.984 / 1000; % mass (kg)
% About the center of mass
Ib = [22503.0 0 0; 0 29230.0 0; 0 0 69400.0] / 100 / 1000; % Moment of inertia (kg cm ^ 2)
Gb = [Ib zeros(3); zeros(3) m * eye(3)];

% Parameters
params.P = P;
params.ls_0 = ls_0;
params.lc_0 = lc_0;
params.ls_max = ls_max;
params.lc_max = lc_max;
params.k = k;
params.C1 = C1;
params.C2 = C2;
params.Gb = Gb;