function [Ke, Fe] = element_stiffness_simplified_mitc4(elem_x, elem_thick, elem_vn, elem_v1, elem_v2, E, nu, surface_load)
% ELEMENT_STIFFNESS_SIMPLIFIED_MITC4 Stiffness matrix and load vector
% for the Simplified MITC4+ shell element.
%
% Inputs:
%   elem_x       : 4 x 3 nodal coordinates [x, y, z]
%   elem_thick   : 4 x 1 node thicknesses (or scalar)
%   elem_vn      : 4 x 3 node director vectors Vn
%   elem_v1      : 4 x 3 node vectors V1
%   elem_v2      : 4 x 3 node vectors V2
%   E, nu        : Elastic parameters
%   surface_load : Surface pressure load or gravity vector [px, py, pz] (optional)

if nargin < 8, surface_load = 0; end
if length(elem_thick) == 1, elem_thick = repmat(elem_thick, 4, 1); end

num_dof = 24; % 4 nodes x 6 DOFs (u, v, w, alpha, beta, gamma)
Ke = zeros(num_dof, num_dof);
Fe = zeros(num_dof, 1);

% Characteristic geometry vectors (Eq. 15)
xi  = [-1,  1,  1, -1];
eta = [-1, -1,  1,  1];

xr = zeros(1, 3); xs = zeros(1, 3); xd = zeros(1, 3);
for i = 1:4
    xr = xr + 0.25 * xi(i) * elem_x(i, :);
    xs = xs + 0.25 * eta(i) * elem_x(i, :);
    xd = xd + 0.25 * xi(i) * eta(i) * elem_x(i, :);
end

% Normal vector to mid-surface plane
n_cross = cross(xr, xs);
norm_n = norm(n_cross);
if norm_n > 1e-12
    n_vec = n_cross / norm_n;
else
    n_vec = [0, 0, 1];
end

% Distortion parameters alpha, beta, d (Eq. 24)
alpha_param = dot(n_vec, cross(xd, xs)) / norm_n;
beta_param  = dot(n_vec, cross(xr, xd)) / norm_n;
d_param     = alpha_param^2 + beta_param^2 - 1;

if abs(d_param) < 1e-12, d_param = -1; end

aA = alpha_param * (alpha_param - 1) / (2 * d_param);
aB = alpha_param * (alpha_param + 1) / (2 * d_param);
aC = beta_param  * (beta_param - 1)  / (2 * d_param);
aD = beta_param  * (beta_param + 1)  / (2 * d_param);
aE = 2 * alpha_param * beta_param    / d_param;

% Matrix C (5 x 5) (Eq. 47)
C_mat = [ ...
    1/2 - aA,  1/2 - aB, -aC,       -aD,       -aE;
   -aA,       -aB,        1/2 - aC,  1/2 - aD, -aE;
    0,         0,         0,         0,         1;
    1/2,      -1/2,       0,         0,         0;
    0,         0,         1/2,      -1/2,       0 ];

% Tying points for membrane strain (A, B, C, D, E)
tying_pts_mem = [
    0,  1;  % A
    0, -1;  % B
    1,  0;  % C
   -1,  0;  % D
    0,  0   % E
];

B_mem_tying = zeros(5, num_dof);

for tp = 1:5
    r_tp = tying_pts_mem(tp, 1);
    s_tp = tying_pts_mem(tp, 2);

    [~, dh_dr, dh_ds] = shape_functions(r_tp, s_tp);

    dx_dr = zeros(1, 3); dx_ds = zeros(1, 3);
    for i = 1:4
        dx_dr = dx_dr + dh_dr(i) * elem_x(i, :);
        dx_ds = dx_ds + dh_ds(i) * elem_x(i, :);
    end

    B_tp_rr = zeros(1, num_dof);
    B_tp_ss = zeros(1, num_dof);
    B_tp_rs = zeros(1, num_dof);

    for i = 1:4
        dof_idx = (i-1)*6 + (1:3);
        B_tp_rr(dof_idx) = dh_dr(i) * dx_dr;
        B_tp_ss(dof_idx) = dh_ds(i) * dx_ds;
        B_tp_rs(dof_idx) = 0.5 * (dh_dr(i) * dx_ds + dh_ds(i) * dx_dr);
    end

    if tp == 1
        B_mem_tying(1, :) = B_tp_rr;
    elseif tp == 2
        B_mem_tying(2, :) = B_tp_rr;
    elseif tp == 3
        B_mem_tying(3, :) = B_tp_ss;
    elseif tp == 4
        B_mem_tying(4, :) = B_tp_ss;
    elseif tp == 5
        B_mem_tying(5, :) = B_tp_rs;
    end
end

% Tying points for transverse shear strain (A, B, C, D)
B_shear_tying = zeros(4, num_dof);
tying_pts_shear = [
    0,  1;  % A (e_rt)
    0, -1;  % B (e_rt)
    1,  0;  % C (e_st)
   -1,  0   % D (e_st)
];

for tp = 1:4
    r_tp = tying_pts_shear(tp, 1);
    s_tp = tying_pts_shear(tp, 2);

    [h_tp, dh_dr_tp, dh_ds_tp] = shape_functions(r_tp, s_tp);

    dx_dr = dh_dr_tp' * elem_x;
    dx_ds = dh_ds_tp' * elem_x;

    Vn_tp = zeros(1, 3);
    for i = 1:4
        Vn_tp = Vn_tp + h_tp(i) * elem_thick(i)/2 * elem_vn(i, :);
    end

    B_rt = zeros(1, num_dof);
    B_st = zeros(1, num_dof);

    for i = 1:4
        u_dof = (i-1)*6 + (1:3);
        rot_dof = (i-1)*6 + (4:5);
        a_i = elem_thick(i);
        v1_i = elem_v1(i, :);
        v2_i = elem_v2(i, :);

        B_rt(u_dof) = 0.5 * dh_dr_tp(i) * Vn_tp;
        B_rt(rot_dof(1)) = 0.5 * h_tp(i) * a_i/2 * dot(dx_dr, -v2_i);
        B_rt(rot_dof(2)) = 0.5 * h_tp(i) * a_i/2 * dot(dx_dr,  v1_i);

        B_st(u_dof) = 0.5 * dh_ds_tp(i) * Vn_tp;
        B_st(rot_dof(1)) = 0.5 * h_tp(i) * a_i/2 * dot(dx_ds, -v2_i);
        B_st(rot_dof(2)) = 0.5 * h_tp(i) * a_i/2 * dot(dx_ds,  v1_i);
    end

    if tp == 1 || tp == 2
        B_shear_tying(tp, :) = B_rt;
    else
        B_shear_tying(tp, :) = B_st;
    end
end

% Local Cartesian frame e_1, e_2, e_3
e3_loc = n_vec;
e1_loc = xr / norm(xr);
e2_loc = cross(e3_loc, e1_loc); e2_loc = e2_loc / norm(e2_loc);

C_plane_stress = E / (1 - nu^2) * [
    1,  nu, 0;
    nu, 1,  0;
    0,  0,  (1-nu)/2
];
G_shear = E / (2 * (1 + nu));
G_trans = 5/6 * G_shear;

D_shell = zeros(5, 5);
D_shell(1:3, 1:3) = C_plane_stress;
D_shell(4, 4)     = G_trans;
D_shell(5, 5)     = G_trans;

[mu_adj, ~] = compute_geometry_parameter(elem_x);

gauss_pts = [-1/sqrt(3), 1/sqrt(3)];
gauss_wts = [1, 1];

if mu_adj > 1e-6
    gauss_pts_adj = mu_adj * gauss_pts;
else
    gauss_pts_adj = gauss_pts;
end

for gi = 1:2
    r_g = gauss_pts_adj(gi); w_r = gauss_wts(gi);
    for gj = 1:2
        s_g = gauss_pts_adj(gj); w_s = gauss_wts(gj);
        for gk = 1:2
            t_g = gauss_pts(gk); w_t = gauss_wts(gk);

            [h_g, dh_dr_g, dh_ds_g] = shape_functions(r_g, s_g);

            g_r = zeros(1, 3); g_s = zeros(1, 3); g_t = zeros(1, 3);
            for i = 1:4
                a_i = elem_thick(i);
                vn_i = elem_vn(i, :);
                g_r = g_r + dh_dr_g(i) * elem_x(i, :) + t_g/2 * a_i * dh_dr_g(i) * vn_i;
                g_s = g_s + dh_ds_g(i) * elem_x(i, :) + t_g/2 * a_i * dh_ds_g(i) * vn_i;
                g_t = g_t + 0.5 * a_i * h_g(i) * vn_i;
            end

            J_mat = [g_r; g_s; g_t];
            detJ = det(J_mat);

            g_hat_r = zeros(1, 3); g_hat_s = zeros(1, 3);
            for i = 1:4
                g_hat_r = g_hat_r + dh_dr_g(i) * elem_x(i, :);
                g_hat_s = g_hat_s + dh_ds_g(i) * elem_x(i, :);
            end

            g_hat_mat = [g_hat_r; g_hat_s; n_vec];
            g_hat_inv = inv(g_hat_mat);
            g_hat_r_contra = g_hat_inv(:, 1)';
            g_hat_s_contra = g_hat_inv(:, 2)';

            alpha_s = 1 + alpha_param * s_g;
            beta_r  = 1 + beta_param  * r_g;

            Q_mat = [ ...
                alpha_s^2,               (beta_param*s_g)^2,      2*beta_param*s_g*alpha_s;
                (alpha_param*r_g)^2,     beta_r^2,                2*alpha_param*r_g*beta_r;
                alpha_param*r_g*alpha_s, beta_param*s_g*beta_r,   alpha_param*beta_param*r_g*s_g + alpha_s*beta_r ...
            ];

            detJ_0 = norm(cross(xr, xs));
            detJ_rs = norm(cross(g_hat_r, g_hat_s));
            lambda_rs = detJ_0 / max(detJ_rs, 1e-12);

            M_mat = [ ...
                1 - 2*alpha_param*lambda_rs*s_g, 0,                               -2*beta_param*lambda_rs*s_g, lambda_rs*s_g, 0;
                0,                               1 - 2*beta_param*lambda_rs*r_g, -2*alpha_param*lambda_rs*r_g, 0,            lambda_rs*r_g;
                0,                               0,                               1,                           0,            0 ...
            ];

            B_mem_assumed = Q_mat * M_mat * C_mat * B_mem_tying;

            B_bend = zeros(3, num_dof);
            for i = 1:4
                rot_dof = (i-1)*6 + (4:5);
                a_i = elem_thick(i);
                v1_i = elem_v1(i, :);
                v2_i = elem_v2(i, :);

                B_bend(1, rot_dof(1)) = t_g * dh_dr_g(i) * a_i/2 * dot(g_hat_r, -v2_i);
                B_bend(1, rot_dof(2)) = t_g * dh_dr_g(i) * a_i/2 * dot(g_hat_r,  v1_i);

                B_bend(2, rot_dof(1)) = t_g * dh_ds_g(i) * a_i/2 * dot(g_hat_s, -v2_i);
                B_bend(2, rot_dof(2)) = t_g * dh_ds_g(i) * a_i/2 * dot(g_hat_s,  v1_i);

                B_bend(3, rot_dof(1)) = 0.5 * t_g * a_i/2 * (dh_dr_g(i)*dot(g_hat_s, -v2_i) + dh_ds_g(i)*dot(g_hat_r, -v2_i));
                B_bend(3, rot_dof(2)) = 0.5 * t_g * a_i/2 * (dh_dr_g(i)*dot(g_hat_s,  v1_i) + dh_ds_g(i)*dot(g_hat_r,  v1_i));
            end

            B_inplane_cov = B_mem_assumed + B_bend;

            B_shear_cov = zeros(2, num_dof);
            B_shear_cov(1, :) = 0.5 * (1 + s_g) * B_shear_tying(1, :) + 0.5 * (1 - s_g) * B_shear_tying(2, :);
            B_shear_cov(2, :) = 0.5 * (1 + r_g) * B_shear_tying(3, :) + 0.5 * (1 - r_g) * B_shear_tying(4, :);

            c11 = dot(g_hat_r_contra, e1_loc); c12 = dot(g_hat_r_contra, e2_loc);
            c21 = dot(g_hat_s_contra, e1_loc); c22 = dot(g_hat_s_contra, e2_loc);

            T_inplane = zeros(3, 3);
            T_inplane(1, 1) = c11^2;      T_inplane(1, 2) = c21^2;      T_inplane(1, 3) = 2*c11*c21;
            T_inplane(2, 1) = c12^2;      T_inplane(2, 2) = c22^2;      T_inplane(2, 3) = 2*c12*c22;
            T_inplane(3, 1) = c11*c12;    T_inplane(3, 2) = c21*c22;    T_inplane(3, 3) = c11*c22 + c12*c21;

            c33 = dot(n_vec, e3_loc);
            T_shear = zeros(2, 2);
            T_shear(1, 1) = c11*c33; T_shear(1, 2) = c21*c33;
            T_shear(2, 1) = c12*c33; T_shear(2, 2) = c22*c33;

            B_local = zeros(5, num_dof);
            B_local(1:3, :) = T_inplane * B_inplane_cov;
            B_local(4:5, :) = 2 * T_shear * B_shear_cov;

            dV = detJ * w_r * w_s * w_t;
            Ke = Ke + B_local' * D_shell * B_local * dV;

            % Project surface/gravity loads into 3D global DOFs using surface normal n_vec
            if gk == 1 && (isnumeric(surface_load) && norm(surface_load) > 0 || isstruct(surface_load))
                p_vec = zeros(1, 3);
                if isnumeric(surface_load)
                    if length(surface_load) == 1
                        p_vec = surface_load * n_vec; % normal pressure
                    elseif length(surface_load) == 3
                        p_vec = surface_load;        % global load vector (e.g. gravity)
                    end
                end

                for i = 1:4
                    u_dof = (i-1)*6 + (1:3);
                    Fe(u_dof) = Fe(u_dof) + (h_g(i) * detJ_rs * w_r * w_s) * p_vec';
                end
            end
        end
    end
end

% Small artificial drilling stiffness for unconstrained 6th DOFs
for i = 1:4
    gamma_dof = (i-1)*6 + 6;
    if Ke(gamma_dof, gamma_dof) < 1e-8
        Ke(gamma_dof, gamma_dof) = 1e-4 * max(diag(Ke));
    end
end

end
