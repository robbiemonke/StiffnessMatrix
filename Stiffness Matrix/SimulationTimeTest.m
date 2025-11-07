clear
clc

% Constants
P = [-3.0162, 0, 0, 1;...
    -1.524, -2.6035, 0, 1;...
    1.524, -2.6035, 0, 1;...
    3.0162, 0, 0, 1]';
ls_0 = 3.175;
lc_0 = 3.810;
k = 0.77;
epsilon = 1e-6;

% Helper Function Definitions
wrench = @(xi) getWrench(xi, P, ls_0, lc_0, k);
TM = @(xi) screw2TMat(xi);

% Minimal energy configuration
T = TMatExponential3([0 1 0]', pi/2, [0 0 0]') * TMatExponential3([0 0 1]', pi, [0 0 0]');
T_translation = [1 0 0 0;...
                 0 1 0 1.57;...
                 0 0 1 0;...
                 0 0 0 1];
T = T * T_translation;
[u, theta, vtilde] = EquivalentScrew3(T);
xi_0 = [u; vtilde] * theta;
wrench_0 = wrench(xi_0);

N = 1000; % number of simulation configurations
% Simulation time and error testing
for i = 1:N
    % Development of random configuration using minimization of energy
    Fext(:, i) = rand(3, 1) * 0.1;
    xi_min(:, i) = fminunc(@(xi) fFindFunction_AA(xi, ls_0, lc_0, k, P, Fext(:, i)), xi_0);

    % Calculation of change in screw and wrench at the configuration
    dxi = xi_min(:, i) - xi_0;
    dwrench = wrench(xi_min(:, i)) - wrench_0;

    % Calculation and timing of numerical stiffness matrix solution
    tic
    for j = 1:6
        dx = zeros(6, 1); 
        dx(j) = epsilon;
        wrench_plus = wrench(xi_min(:, i) + dx);
        wrench_minus = wrench(xi_min(:, i) - dx);
        K_fd(:, j, i) = (wrench_plus - wrench_minus) / (2 * epsilon);
    end
    num_SM_time(i) = toc;
    dwrench_num(:, :, i) = K_fd(:, :, i) * dxi;

    % Calculation and timing of analytical stiffness matrix solution
    tic
    K_ana(:, :, i) = getStiffnessMatrix(xi_min(:, i), P, ls_0, lc_0, k);
    ana_SM_time(i) = toc;
    dwrench_ana(:, :, i) = K_ana(:, :, i) * dxi;

    % Calculation and timing of numerical transformation matrix derivative solution
    tic
    for j = 1:6
        dx = zeros(6, 1); 
        dx(j) = epsilon;
        T_plus = TM(xi_min(:, i) + dx);
        T_minus = TM(xi_min(:, i) - dx);
        dT_num(:, :, j) = (T_plus - T_minus) / (2 * epsilon);
    end
    num_TM_time(i) = toc;


    % Calculation and timing of analytical transformation matrix derivative solution
    tic
    dT_ana = screw2dT(xi_min(:, i));
    ana_TM_time(i) = toc;
end

% Display of results
disp("Average analytical time to solve 1000 random tests SM (ms):")
disp(mean(ana_SM_time) * 1000)
disp("Average numerical time to solve 1000 random tests SM (ms):")
disp(mean(num_SM_time) * 1000)
disp("Average error SM:")
disp(norm(dwrench_ana - dwrench_num, "fro") / norm(dwrench_ana, "fro"))

disp("Average analytical time to solve 1000 random tests TM (ms):")
disp(mean(ana_TM_time) * 1000)
disp("Average numerical time to solve 1000 random tests TM (ms):")
disp(mean(num_TM_time) * 1000)
disp("Average error TM:")
disp(norm(dT_ana - dT_num, "fro") / norm(dT_ana, "fro"))