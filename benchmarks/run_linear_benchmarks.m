function run_linear_benchmarks()
% RUN_LINEAR_BENCHMARKS Solves linear benchmark problems from Section 5:
% - L-shaped structure
% - Scordelis-Lo roof

addpath('../src'); addpath('../utils');

fprintf('====================================================\n');
fprintf('  RUNNING LINEAR BENCHMARK PROBLEMS (SECTION 5)\n');
fprintf('====================================================\n\n');

%% 1. L-SHAPED STRUCTURE (Section 5.1)
fprintf('[1] L-Shaped Structure Problem:\n');
L = 50.0; E = 2.0e11; nu = 0.22; p = 1.0; thickness = 0.6;
mesh_N = [4, 8, 16];
errors_reg = zeros(1, length(mesh_N));

% Compute reference strain energy with fine L-shaped mesh N = 32
[ref_nodes, ref_elems] = generate_lshape_mesh(L, 32);
bc_ref = [];
% Fixed bottom edge (y = 0)
for i = 1:size(ref_nodes, 1)
    if abs(ref_nodes(i, 2)) < 1e-5
        for dof = 1:6, bc_ref = [bc_ref; i, dof, 0]; end
    end
end

[U_ref, ~, K_ref, ~] = solve_linear_fe(ref_nodes, ref_elems, thickness, E, nu, bc_ref, struct('surface_pressure', p));
energy_ref = 0.5 * U_ref' * K_ref * U_ref;

for idx = 1:length(mesh_N)
    N = mesh_N(idx);
    [nodes, elements] = generate_lshape_mesh(L, N);

    bc = [];
    for i = 1:size(nodes, 1)
        if abs(nodes(i, 2)) < 1e-5
            for dof = 1:6, bc = [bc; i, dof, 0]; end
        end
    end

    [U, ~, K, ~] = solve_linear_fe(nodes, elements, thickness, E, nu, bc, struct('surface_pressure', p));
    energy_h = 0.5 * U' * K * U;

    % Relative strain energy norm error Eh
    errors_reg(idx) = sqrt(abs(energy_ref - energy_h) / energy_ref);
    fprintf('    Mesh N=%2d: Relative Energy Error Eh = %.6e\n', N, errors_reg(idx));
end

plot_convergence(mesh_N, errors_reg, {'Simplified MITC4+'}, 'L-Shaped Structure Energy Error Convergence');

%% 2. SCORDELIS-LO ROOF (Section 5.3)
fprintf('\n[2] Scordelis-Lo Roof Problem:\n');
% Scordelis-Lo roof geometry: L = 50.0, R = 25.0, theta = 40 deg = 40*pi/180
L = 50.0; R = 25.0; phi_span = 40 * pi / 180; E = 4.32e8; nu = 0.0; thickness = 0.25; fz = -90.0;
[nodes_roof, elements_roof] = generate_quad_mesh([0, L], [-phi_span, phi_span], 8, 8, 'none');

% Convert cylindrical coords (r=R, theta) to Cartesian
for i = 1:size(nodes_roof, 1)
    y_val = nodes_roof(i, 1);
    theta = nodes_roof(i, 2);
    nodes_roof(i, :) = [R * sin(theta), y_val, R * cos(theta)];
end

% Fixed rigid diaphragm ends at y = 0 and y = L (u_x = u_z = 0)
% Fixed axial displacement u_y = 0 at y = L/2 symmetry line
bc_roof = [];
for i = 1:size(nodes_roof, 1)
    y_i = nodes_roof(i, 2);
    if abs(y_i) < 1e-5 || abs(y_i - L) < 1e-5
        bc_roof = [bc_roof; i, 1, 0; i, 3, 0];
    end
    if abs(y_i - L/2) < 1e-3
        bc_roof = [bc_roof; i, 2, 0];
    end
end

% Self-weight gravity load in Z direction per unit surface area passed in struct
loads_roof = struct('surface_load', [0, 0, fz]);
[U_roof, ~, ~, ~] = solve_linear_fe(nodes_roof, elements_roof, thickness, E, nu, bc_roof, loads_roof);

% Find mid-point node at x=0, y=L/2
mid_node = 1; min_d = 1e9;
for i = 1:size(nodes_roof, 1)
    d = norm(nodes_roof(i, :) - [0, L/2, R]);
    if d < min_d
        min_d = d; mid_node = i;
    end
end

w_mid = U_roof((mid_node-1)*6 + 3);
fprintf('    Mid-point vertical displacement w_C: %.4f (Analytical Reference: -0.3024)\n', w_mid);
plot_deformation(nodes_roof, elements_roof, U_roof, 5.0, 'Scordelis-Lo Roof Deformation');

end
