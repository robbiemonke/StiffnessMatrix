clear; 
clc;
close all;

%% Constants
constants;

%% Optimization
% Initial guess for transformation matrix and screw
T = TMatExponential3([0 1 0]', pi/2, [0 0 0]') * TMatExponential3([0 0 1]', pi, [0 0 0]');
T_translation = [1 0 0 0;...
                 0 1 0 2;...
                 0 0 1 0;...
                 0 0 0 1];
T = T * T_translation;
[u, theta, vtilde] = EquivalentScrew3(T);
xi_0 = [u; vtilde] * theta;

Fext = 0 * [0 0 -1]';

% Optimization of screw using minimization of energy functions. 
options = optimoptions(@fminunc,'Display','iter');
xi_min = fminunc(@(xi) fFindFunction(xi, params, Fext), xi_0, options);
T_min = screw2TMat(xi_min);

%% Wrench
wrench = getWrench(xi_min, params);

%% Plotting model as semi-circles
% Development of points for plotting
P1 = P;
P2 = T_min * P;

A1 = P1(1:3, 1);
B1 = P1(1:3, 2);
C1 = P1(1:3, 3);
D1 = P1(1:3, 4);

A2 = P2(1:3, 1);
B2 = P2(1:3, 2);
C2 = P2(1:3, 3);
D2 = P2(1:3, 4);

% Finding Euler angles from rotation matrix
eulerZYX = rotm2eul(T_min(1:3, 1:3), "XYZ");

% Plotting of arc points - nodes and semicircle arc
figure(1)
scatter3(P1(1, :), P1(2, :), P1(3, :))
text(P1(1, :), P1(2, :), P1(3, :)+0.5, ["A1", "B1", "C1", "D1"], "Color", "r")
hold on
[semi1] = semiCirclePoints3D(P(1,4), 101, [0, 0, 0], 0, 0, 0);
plot3(semi1(:,1), semi1(:,2), semi1(:,3), "Color", "b", "LineWidth", 2)

scatter3(P2(1, :), P2(2, :), P2(3, :))
text(P2(1, :), P2(2, :), P2(3, :)+0.5, ["A2", "B2", "C2", "D2"], "Color", "r")
[semi2] = semiCirclePoints3D(P(1,4), 101, T_min(1:3, 4)', eulerZYX(1), eulerZYX(2), eulerZYX(3));
plot3(semi2(:,1), semi2(:,2), semi2(:,3), "Color", "r", "LineWidth", 2)

% Plot of diagonal strings.
plot3([A1(1), A2(1)], [A1(2), A2(2)], [A1(3), A2(3)])
plot3([A1(1), D2(1)], [A1(2), D2(2)], [A1(3), D2(3)])
plot3([D1(1), A2(1)], [D1(2), A2(2)], [D1(3), A2(3)])
plot3([D1(1), D2(1)], [D1(2), D2(2)], [D1(3), D2(3)])

% Plot of straight strings
plot3([A1(1), B2(1)], [A1(2), B2(2)], [A1(3), B2(3)])
plot3([A1(1), C2(1)], [A1(2), C2(2)], [A1(3), C2(3)])
plot3([B1(1), A2(1)], [B1(2), A2(2)], [B1(3), A2(3)])
plot3([B1(1), D2(1)], [B1(2), D2(2)], [B1(3), D2(3)])
plot3([C1(1), A2(1)], [C1(2), A2(2)], [C1(3), A2(3)])
plot3([C1(1), D2(1)], [C1(2), D2(2)], [C1(3), D2(3)])
plot3([D1(1), B2(1)], [D1(2), B2(2)], [D1(3), B2(3)])
plot3([D1(1), C2(1)], [D1(2), C2(2)], [D1(3), C2(3)])

%% Plotting model from STL file
figure(2)

% Imports stl
arc1 = fegeometry("arc1-CenterHole.stl");

% Scales from mm to cm
arc1 = scale(arc1, 1/10);

% Rotates bodies to align with points
arc1 = rotate(arc1, -90, [0, 0, 0], [0 0 1]);
arc1 = rotate(arc1, 90, [0 0 0], [0 1 0]);
arc2 = rotate(arc1, eulerZYX(1) * (180/pi), [0 0 0], [1 0 0]);
arc2 = rotate(arc2, eulerZYX(2) * (180/pi), [0 0 0], [0 1 0]);
arc2 = rotate(arc2, eulerZYX(3) * (180/pi), [0 0 0], [0 0 1]);
arc2 = translate(arc2, T_min(1:3, 4)');

% Plots both arcs
pdegplot(arc1)
hold on
pdegplot(arc2)
delete(findobj(gca,'type','Text')); 
delete(findobj(gca,'type','Quiver')); 
hold on

text(P1(1, :), P1(2, :), P1(3, :) + 1, ["A1", "B1", "C1", "D1"], "Color", "r")
text(P2(1, :) + 1, P2(2, :), P2(3, :), ["A2", "B2", "C2", "D2"], "Color", "r")

% Plot of diagonal strings.
plot3([A1(1), A2(1)], [A1(2), A2(2)], [A1(3), A2(3)])
plot3([A1(1), D2(1)], [A1(2), D2(2)], [A1(3), D2(3)])
plot3([D1(1), A2(1)], [D1(2), A2(2)], [D1(3), A2(3)])
plot3([D1(1), D2(1)], [D1(2), D2(2)], [D1(3), D2(3)])

% Plot of straight strings
plot3([A1(1), B2(1)], [A1(2), B2(2)], [A1(3), B2(3)])
plot3([A1(1), C2(1)], [A1(2), C2(2)], [A1(3), C2(3)])
plot3([B1(1), A2(1)], [B1(2), A2(2)], [B1(3), A2(3)])
plot3([B1(1), D2(1)], [B1(2), D2(2)], [B1(3), D2(3)])
plot3([C1(1), A2(1)], [C1(2), A2(2)], [C1(3), A2(3)])
plot3([C1(1), D2(1)], [C1(2), D2(2)], [C1(3), D2(3)])
plot3([D1(1), B2(1)], [D1(2), B2(2)], [D1(3), B2(3)])
plot3([D1(1), C2(1)], [D1(2), C2(2)], [D1(3), C2(3)])

%% Plotting of force applied to model
% quiver3(0, 6, 0, 0, -2.5, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 3)
% text(0, 6, 0, "Perturbation Force", "Color", "black", 'FontSize', 28)

% quiver3(A2(1)+0.32, A2(2), A2(3), 0, 0, -2.5, 'r', 'LineWidth', 3, 'MaxHeadSize', 3)
% text(A2(1), A2(2), A2(3)-3, "Perturbation Force", "Color", "black", 'FontSize', 28)