function [U_history, load_factors, U_final] = solve_nonlinear_fe(nodes, elements, thickness, E, nu, boundary_conditions, max_load, num_steps)
% SOLVE_NONLINEAR_FE Incremental-Iterative Newton-Raphson solver for
% geometric nonlinear analysis using Simplified MITC4+ shell elements.

if nargin < 8, num_steps = 10; end

num_nodes = size(nodes, 1);
total_dofs = num_nodes * 6;
num_elem  = size(elements, 1);

[Vn, V1, V2] = compute_director_vectors(nodes, elements);

U_current = zeros(total_dofs, 1);
U_history = zeros(total_dofs, num_steps + 1);
load_factors = linspace(0, 1, num_steps + 1);

% Fixed DOFs
fixed_dofs = [];
if ~isempty(boundary_conditions)
    if size(boundary_conditions, 2) == 2
        fixed_dofs = boundary_conditions(:, 1);
    elseif size(boundary_conditions, 2) == 3
        for b = 1:size(boundary_conditions, 1)
            dof_idx = (boundary_conditions(b, 1) - 1) * 6 + boundary_conditions(b, 2);
            fixed_dofs = [fixed_dofs; dof_idx];
        end
    end
end
free_dofs = setdiff(1:total_dofs, fixed_dofs);

% Target total external load vector at full load (lambda = 1.0)
F_ext_full = zeros(total_dofs, 1);
if isstruct(max_load) && isfield(max_load, 'point_loads')
    p_loads = max_load.point_loads;
    for p = 1:size(p_loads, 1)
        dof_idx = (p_loads(p, 1) - 1) * 6 + p_loads(p, 2);
        F_ext_full(dof_idx) = F_ext_full(dof_idx) + p_loads(p, 3);
    end
end

% Step-by-step loading
for step = 1:num_steps
    lambda = load_factors(step + 1);
    F_ext_step = lambda * F_ext_full;

    % Newton-Raphson iterations
    max_iter = 25;
    tol = 1e-4;

    for iter = 1:max_iter
        % Assemble global tangent stiffness and internal force
        K_tangent = zeros(total_dofs, total_dofs);
        F_internal = zeros(total_dofs, 1);

        for e = 1:num_elem
            elem_nodes = elements(e, :);
            elem_x = nodes(elem_nodes, :);
            if length(thickness) == 1, elem_t = thickness; else, elem_t = thickness(elem_nodes); end

            elem_dofs = zeros(1, 24);
            for i = 1:4, elem_dofs((i-1)*6 + (1:6)) = (elem_nodes(i)-1)*6 + (1:6); end
            elem_u = U_current(elem_dofs);

            [Ke_nl, Fe_int] = element_nonlinear_simplified_mitc4(elem_x, elem_u, elem_t, Vn(elem_nodes,:), V1(elem_nodes,:), V2(elem_nodes,:), E, nu);

            K_tangent(elem_dofs, elem_dofs) = K_tangent(elem_dofs, elem_dofs) + Ke_nl;
            F_internal(elem_dofs) = F_internal(elem_dofs) + Fe_int;
        end

        % Residual vector
        R = F_ext_step - F_internal;
        R(fixed_dofs) = 0;

        norm_R = norm(R(free_dofs));
        if norm_R < tol * max(norm(F_ext_step(free_dofs)), 1.0)
            break;
        end

        % Displacement increment
        dU = zeros(total_dofs, 1);
        dU(free_dofs) = K_tangent(free_dofs, free_dofs) \ R(free_dofs);

        U_current = U_current + dU;
    end

    U_history(:, step + 1) = U_current;
end

U_final = U_current;

end
