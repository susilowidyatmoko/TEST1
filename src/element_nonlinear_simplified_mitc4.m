function [Ke_nl, Finternal] = element_nonlinear_simplified_mitc4(elem_x, elem_u, elem_thick, elem_vn, elem_v1, elem_v2, E, nu)
% ELEMENT_NONLINEAR_SIMPLIFIED_MITC4 Total Lagrangian geometric nonlinear
% formulation for Simplified MITC4+ shell element.
%
% Calculates Green-Lagrange strain, 2nd Piola-Kirchhoff stress,
% internal force vector Finternal, and tangent stiffness matrix Ke_nl.

if length(elem_thick) == 1, elem_thick = repmat(elem_thick, 4, 1); end

num_dof = 24;
Ke_nl = zeros(num_dof, num_dof);
Finternal = zeros(num_dof, 1);

% Linear stiffness & strain-displacement relation matrix
[Ke_lin, ~] = element_stiffness_simplified_mitc4(elem_x, elem_thick, elem_vn, elem_v1, elem_v2, E, nu, 0);

% Material law matrix D
C_plane_stress = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
G_trans = 5/6 * (E / (2 * (1 + nu)));
D_shell = zeros(5, 5);
D_shell(1:3, 1:3) = C_plane_stress;
D_shell(4, 4) = G_trans;
D_shell(5, 5) = G_trans;

% Integration points
gauss_pts = [-1/sqrt(3), 1/sqrt(3)];
gauss_wts = [1, 1];
[mu_adj, ~] = compute_geometry_parameter(elem_x);
if mu_adj > 1e-6, gauss_pts_adj = mu_adj * gauss_pts; else, gauss_pts_adj = gauss_pts; end

% Mid-surface geometry
xi  = [-1,  1,  1, -1];
eta = [-1, -1,  1,  1];
xr = zeros(1, 3); xs = zeros(1, 3);
for i = 1:4
    xr = xr + 0.25 * xi(i) * elem_x(i, :);
    xs = xs + 0.25 * eta(i) * elem_x(i, :);
end
n_vec = cross(xr, xs); n_vec = n_vec / norm(n_vec);
e1_loc = xr / norm(xr);
e2_loc = cross(n_vec, e1_loc); e2_loc = e2_loc / norm(e2_loc);

Kg = zeros(num_dof, num_dof);

for gi = 1:2
    r_g = gauss_pts_adj(gi); w_r = gauss_wts(gi);
    for gj = 1:2
        s_g = gauss_pts_adj(gj); w_s = gauss_wts(gj);
        for gk = 1:2
            t_g = gauss_pts(gk); w_t = gauss_wts(gk);

            [h_g, dh_dr_g, dh_ds_g] = shape_functions(r_g, s_g);

            [detJ, ~, B_local] = evaluate_element_B(elem_x, elem_thick, elem_vn, elem_v1, elem_v2, r_g, s_g, t_g, xr, xs, n_vec, e1_loc, e2_loc);
            dV = detJ * w_r * w_s * w_t;

            % Compute 3D displacement gradient matrix H = du_i / dx_j (3 x 3)
            dx_dr = zeros(1, 3); dx_ds = zeros(1, 3);
            for i = 1:4
                dx_dr = dx_dr + dh_dr_g(i) * elem_x(i, :);
                dx_ds = dx_ds + dh_ds_g(i) * elem_x(i, :);
            end

            u_nodes = zeros(4, 3);
            for i = 1:4
                u_nodes(i, :) = elem_u((i-1)*6 + (1:3))';
            end

            du_dr = dh_dr_g' * u_nodes;
            du_ds = dh_ds_g' * u_nodes;

            H_local = zeros(3, 3);
            H_local(:, 1) = du_dr';
            H_local(:, 2) = du_ds';

            e_linear = B_local * elem_u;

            % 2nd Piola-Kirchhoff stress S
            sigma = D_shell * e_linear;

            % Internal force accumulation: Finternal = integral (B_L' * S dV)
            Finternal = Finternal + B_local' * sigma * dV;

            % Geometric stiffness matrix Kg = integral (G' * S_block * G dV)
            G_mat = zeros(9, num_dof);
            for i = 1:4
                u_dof = (i-1)*6 + (1:3);
                G_mat(1, u_dof(1)) = dh_dr_g(i);
                G_mat(2, u_dof(2)) = dh_dr_g(i);
                G_mat(3, u_dof(3)) = dh_dr_g(i);

                G_mat(4, u_dof(1)) = dh_ds_g(i);
                G_mat(5, u_dof(2)) = dh_ds_g(i);
                G_mat(6, u_dof(3)) = dh_ds_g(i);

                G_mat(7, u_dof(1)) = h_g(i);
                G_mat(8, u_dof(2)) = h_g(i);
                G_mat(9, u_dof(3)) = h_g(i);
            end

            S_3x3 = [
                sigma(1), sigma(3), 0;
                sigma(3), sigma(2), 0;
                0,        0,        sigma(4)
            ];

            S_block = kron(S_3x3, eye(3));
            Kg = Kg + G_mat' * S_block * G_mat * dV;
        end
    end
end

Ke_nl = Ke_lin + Kg;

end

function [detJ, dV, B_local] = evaluate_element_B(elem_x, elem_thick, elem_vn, elem_v1, elem_v2, r_g, s_g, t_g, xr, xs, n_vec, e1_loc, e2_loc)
num_dof = 24;
[h_g, dh_dr_g, dh_ds_g] = shape_functions(r_g, s_g);

g_r = zeros(1, 3); g_s = zeros(1, 3); g_t = zeros(1, 3);
for i = 1:4
    a_i = elem_thick(i); vn_i = elem_vn(i, :);
    g_r = g_r + dh_dr_g(i) * elem_x(i, :) + t_g/2 * a_i * dh_dr_g(i) * vn_i;
    g_s = g_s + dh_ds_g(i) * elem_x(i, :) + t_g/2 * a_i * dh_ds_g(i) * vn_i;
    g_t = g_t + 0.5 * a_i * h_g(i) * vn_i;
end

detJ = det([g_r; g_s; g_t]);
dV = detJ;

g_hat_r = zeros(1, 3); g_hat_s = zeros(1, 3);
for i = 1:4
    g_hat_r = g_hat_r + dh_dr_g(i) * elem_x(i, :);
    g_hat_s = g_hat_s + dh_ds_g(i) * elem_x(i, :);
end
g_hat_mat = [g_hat_r; g_hat_s; n_vec];
g_hat_inv = inv(g_hat_mat);
g_hat_r_contra = g_hat_inv(:, 1)';
g_hat_s_contra = g_hat_inv(:, 2)';

c11 = dot(g_hat_r_contra, e1_loc); c12 = dot(g_hat_r_contra, e2_loc);
c21 = dot(g_hat_s_contra, e1_loc); c22 = dot(g_hat_s_contra, e2_loc);

T_inplane = zeros(3, 3);
T_inplane(1, 1) = c11^2;      T_inplane(1, 2) = c21^2;      T_inplane(1, 3) = 2*c11*c21;
T_inplane(2, 1) = c12^2;      T_inplane(2, 2) = c22^2;      T_inplane(2, 3) = 2*c12*c22;
T_inplane(3, 1) = c11*c12;    T_inplane(3, 2) = c21*c22;    T_inplane(3, 3) = c11*c22 + c12*c21;

c33 = dot(n_vec, n_vec);
T_shear = [c11*c33, c21*c33; c12*c33, c22*c33];

B_inplane_cov = zeros(3, num_dof);
B_shear_cov = zeros(2, num_dof);

for i = 1:4
    u_dof = (i-1)*6 + (1:3);
    rot_dof = (i-1)*6 + (4:5);
    a_i = elem_thick(i);
    v1_i = elem_v1(i, :); v2_i = elem_v2(i, :);

    B_inplane_cov(1, u_dof) = dh_dr_g(i) * g_hat_r;
    B_inplane_cov(2, u_dof) = dh_ds_g(i) * g_hat_s;
    B_inplane_cov(3, u_dof) = 0.5 * (dh_dr_g(i)*g_hat_s + dh_ds_g(i)*g_hat_r);

    B_inplane_cov(1, rot_dof(1)) = t_g * dh_dr_g(i) * a_i/2 * dot(g_hat_r, -v2_i);
    B_inplane_cov(1, rot_dof(2)) = t_g * dh_dr_g(i) * a_i/2 * dot(g_hat_r,  v1_i);
    B_inplane_cov(2, rot_dof(1)) = t_g * dh_ds_g(i) * a_i/2 * dot(g_hat_s, -v2_i);
    B_inplane_cov(2, rot_dof(2)) = t_g * dh_ds_g(i) * a_i/2 * dot(g_hat_s,  v1_i);

    B_shear_cov(1, u_dof) = 0.5 * dh_dr_g(i) * n_vec;
    B_shear_cov(2, u_dof) = 0.5 * dh_ds_g(i) * n_vec;
end

B_local = zeros(5, num_dof);
B_local(1:3, :) = T_inplane * B_inplane_cov;
B_local(4:5, :) = 2 * T_shear * B_shear_cov;

end
