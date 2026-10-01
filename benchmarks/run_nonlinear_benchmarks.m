function run_nonlinear_benchmarks()
% RUN_NONLINEAR_BENCHMARKS Solves geometric nonlinear benchmark problems:
% - Large deflection cantilever beam/plate
% - Hemispherical shell with hole (Section 5.7)

addpath('../src'); addpath('../utils');

fprintf('====================================================\n');
fprintf('  RUNNING GEOMETRIC NONLINEAR BENCHMARKS (SECTION 5)\n');
fprintf('====================================================\n\n');

%% 1. CANTILEVER BEAM UNDER TIP LOAD
fprintf('[1] Large Deflection Cantilever Beam/Plate:\n');
L = 10.0; W = 1.0; thickness = 0.1; E = 1.2e6; nu = 0.0;
[nodes, elements] = generate_quad_mesh([0, L], [0, W], 10, 2, 'none');

% Clamped at x = 0
bc = [];
for i = 1:size(nodes, 1)
    if abs(nodes(i, 1)) < 1e-6
        for dof = 1:6, bc = [bc; i, dof, 0]; end
    end
end

% Tip vertical load P at x = L
tip_nodes = find(abs(nodes(:, 1) - L) < 1e-6);
P_total = 4.0;
P_node = P_total / length(tip_nodes);

p_loads = [];
for t = 1:length(tip_nodes)
    p_loads = [p_loads; tip_nodes(t), 3, P_node];
end
max_load = struct('point_loads', p_loads);

num_steps = 10;
[U_hist, load_factors, U_final] = solve_nonlinear_fe(nodes, elements, thickness, E, nu, bc, max_load, num_steps);

w_tip = zeros(1, num_steps + 1);
for s = 1:num_steps + 1
    w_tip(s) = U_hist((tip_nodes(1)-1)*6 + 3, s);
end

fprintf('    Final tip vertical displacement w_tip: %.4f\n', w_tip(end));

figure('Color', 'w', 'Name', 'Nonlinear Load-Displacement Curve');
plot(w_tip, load_factors * P_total, 'b-o', 'LineWidth', 1.8, 'MarkerFaceColor', 'b');
grid on; box on;
xlabel('Tip Vertical Displacement w'); ylabel('Applied Load P');
title('Large Deflection Cantilever Load-Displacement', 'FontSize', 12, 'FontWeight', 'bold');

plot_deformation(nodes, elements, U_final, 1.0, 'Nonlinear Deformed Cantilever');

%% 2. HEMISPHERICAL SHELL WITH HOLE (Section 5.7)
fprintf('\n[2] Hemispherical Shell with Hole (1/4 Symmetry):\n');
R = 10.0; phi_0 = 18 * pi/180; thickness = 0.04; E = 6.825e7; nu = 0.3; P_max = 400.0;
N_phi = 8; N_theta = 8;

[nodes_hemi, elements_hemi] = generate_quad_mesh([phi_0, pi/2], [0, pi/2], N_phi, N_theta, 'none');

% Spherical coordinates mapping
for i = 1:size(nodes_hemi, 1)
    phi   = nodes_hemi(i, 1);
    theta = nodes_hemi(i, 2);
    nodes_hemi(i, :) = [R * sin(phi) * cos(theta), R * sin(phi) * sin(theta), R * cos(phi)];
end

% Symmetry BCs on x-z plane (theta = 0) and y-z plane (theta = pi/2)
bc_hemi = [];
for i = 1:size(nodes_hemi, 1)
    if abs(nodes_hemi(i, 2)) < 1e-5 % x-z plane: u_y = 0, rot_x = 0, rot_z = 0
        bc_hemi = [bc_hemi; i, 2, 0; i, 4, 0; i, 6, 0];
    end
    if abs(nodes_hemi(i, 1)) < 1e-5 % y-z plane: u_x = 0, rot_y = 0, rot_z = 0
        bc_hemi = [bc_hemi; i, 1, 0; i, 5, 0; i, 6, 0];
    end
end

% Constrain Z-rigid body motion at apex node (phi = phi_0, top rim)
apex_node = 1; min_d = 1e9;
for i = 1:size(nodes_hemi, 1)
    d = norm(nodes_hemi(i, :) - [0, 0, R*cos(phi_0)]);
    if d < min_d
        min_d = d; apex_node = i;
    end
end
bc_hemi = [bc_hemi; apex_node, 3, 0];

% Point load at A (theta = 0, phi = pi/2) and B (theta = pi/2, phi = pi/2)
node_A = find(abs(nodes_hemi(:,2)) < 1e-5 & abs(nodes_hemi(:,3)) < 1e-2, 1);
node_B = find(abs(nodes_hemi(:,1)) < 1e-5 & abs(nodes_hemi(:,3)) < 1e-2, 1);

p_loads_hemi = [
    node_A, 1,  P_max; % Pulling force along X
    node_B, 2, -P_max  % Pushing force along Y
];
load_hemi = struct('point_loads', p_loads_hemi);

[U_hist_h, ~, U_final_h] = solve_nonlinear_fe(nodes_hemi, elements_hemi, thickness, E, nu, bc_hemi, load_hemi, 5);

u_A = U_final_h((node_A-1)*6 + 1);
v_B = U_final_h((node_B-1)*6 + 2);

fprintf('    Hemispherical Shell Displacements at P = %.1f:\n', P_max);
fprintf('    Radial displacement at A (u_A): %.4f\n', u_A);
fprintf('    Radial displacement at B (v_B): %.4f\n', v_B);

plot_deformation(nodes_hemi, elements_hemi, U_final_h, 1.0, 'Hemispherical Shell Deformed Geometry');

end
