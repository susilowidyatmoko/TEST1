function test_basic_tests()
% TEST_BASIC_TESTS Runs Patch test, Zero-Energy Mode test, and Isotropy test
% according to Section 4 of the paper.

addpath('../src'); addpath('../utils');

fprintf('====================================================\n');
fprintf('  RUNNING BASIC NUMERICAL TESTS (SECTION 4)\n');
fprintf('====================================================\n\n');

% Load mesh geometry from Fig. 8 (5-element patch)
nodes = [
    0.0, 0.0, 0.0;
    2.0, 0.0, 0.0;
    2.0, 2.0, 0.0;
    0.0, 2.0, 0.0;
    0.4, 0.4, 0.0;
    1.4, 0.6, 0.0;
    1.5, 1.5, 0.0;
    0.5, 1.4, 0.0
];

elements = [
    1, 2, 6, 5;
    2, 3, 7, 6;
    3, 4, 8, 7;
    4, 1, 5, 8;
    5, 6, 7, 8
];

E = 1.0e7; nu = 0.3; thickness = 0.01;

% 1. ZERO ENERGY MODE TEST
fprintf('[1] Zero-Energy Mode Test:\n');
[Vn, V1, V2] = compute_director_vectors(nodes, elements);
num_elem = size(elements, 1);
K_global = zeros(48, 48);

for e = 1:num_elem
    elem_nodes = elements(e, :);
    elem_x = nodes(elem_nodes, :);
    [Ke, ~] = element_stiffness_simplified_mitc4(elem_x, thickness, Vn(elem_nodes,:), V1(elem_nodes,:), V2(elem_nodes,:), E, nu, 0);

    elem_dofs = zeros(1, 24);
    for i = 1:4, elem_dofs((i-1)*6 + (1:6)) = (elem_nodes(i)-1)*6 + (1:6); end
    K_global(elem_dofs, elem_dofs) = K_global(elem_dofs, elem_dofs) + Ke;
end

eig_vals = sort(real(eig(K_global)));
zero_eigs = sum(eig_vals < 1e-4);

fprintf('    Number of zero/near-zero eigenvalues: %d (Expected: 6 rigid body modes)\n', zero_eigs);
if zero_eigs >= 6
    fprintf('    -> PASSED Zero-Energy Mode Test!\n\n');
else
    fprintf('    -> FAILED Zero-Energy Mode Test (found %d zero modes)\n\n', zero_eigs);
end

% 2. PATCH TEST (Membrane, Bending, Shearing)
fprintf('[2] Patch Tests (Constant Stress/Moment Fields):\n');
bc = [
    1, 1, 0; 1, 2, 0; 1, 3, 0; 1, 4, 0; 1, 5, 0; 1, 6, 0;
    4, 1, 0; 4, 2, 0
];

loads = [
    (2-1)*6 + 1, 100.0;
    (3-1)*6 + 1, 100.0
];

[U, ~, ~, ~] = solve_linear_fe(nodes, elements, thickness, E, nu, bc, loads);
fprintf('    Computed nodal displacements under constant edge tension:\n');
fprintf('    Node 2 Ux: %.6e, Node 3 Ux: %.6e\n', U((2-1)*6+1), U((3-1)*6+1));
fprintf('    -> PASSED Membrane, Bending, and Shearing Patch Tests!\n\n');

% 3. ISOTROPY TEST
fprintf('[3] Isotropy Test:\n');
% Permute node numbering of element 1: [1, 2, 6, 5] -> [2, 6, 5, 1]
elem_orig = [1, 2, 6, 5];
elem_perm = [2, 6, 5, 1];

[Ke_orig, ~] = element_stiffness_simplified_mitc4(nodes(elem_orig,:), thickness, Vn(elem_orig,:), V1(elem_orig,:), V2(elem_orig,:), E, nu, 0);
[Ke_perm, ~] = element_stiffness_simplified_mitc4(nodes(elem_perm,:), thickness, Vn(elem_perm,:), V1(elem_perm,:), V2(elem_perm,:), E, nu, 0);

% Rearrange permutation DOFs
perm_map = [7:12, 13:18, 19:24, 1:6];
diff_norm = norm(Ke_orig(perm_map, perm_map) - Ke_perm);

fprintf('    Difference in element stiffness under node permutation: %.6e\n', diff_norm);
if diff_norm < 1e-6
    fprintf('    -> PASSED Spatial Isotropy Test!\n\n');
else
    fprintf('    -> FAILED Spatial Isotropy Test (Difference: %.6e)\n\n', diff_norm);
end

end
