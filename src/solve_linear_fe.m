function [U, R_net, K_global, F_ext] = solve_linear_fe(nodes, elements, thickness, E, nu, boundary_conditions, loads)
% SOLVE_LINEAR_FE Assembles global stiffness matrix and solves linear FE system
%
% Inputs:
%   nodes              : N x 3 nodal coordinates
%   elements           : M x 4 element connectivity matrix
%   thickness          : Shell thickness (scalar or N x 1 vector)
%   E, nu              : Material properties
%   boundary_conditions: K x 2 or K x 3 matrix [node_id, dof_index, (value)]
%   loads              : Scalar (surface load), struct ('surface_pressure'), or matrix [dof_global_index, force_value]

num_nodes = size(nodes, 1);
num_elem  = size(elements, 1);
total_dofs = num_nodes * 6;

% Director vectors
[Vn, V1, V2] = compute_director_vectors(nodes, elements);

% Global stiffness and force
K_global = zeros(total_dofs, total_dofs);
F_ext    = zeros(total_dofs, 1);

% Parse surface load / pressure
surface_pressure = 0;
if isstruct(loads) && isfield(loads, 'surface_pressure')
    surface_pressure = loads.surface_pressure;
elseif isnumeric(loads) && numel(loads) == 1
    surface_pressure = loads;
end

% Assemble global stiffness matrix and element surface loads
for e = 1:num_elem
    elem_nodes = elements(e, :);
    elem_x = nodes(elem_nodes, :);

    if length(thickness) == 1
        elem_t = thickness;
    else
        elem_t = thickness(elem_nodes);
    end

    elem_vn = Vn(elem_nodes, :);
    elem_v1 = V1(elem_nodes, :);
    elem_v2 = V2(elem_nodes, :);

    [Ke, Fe] = element_stiffness_simplified_mitc4(elem_x, elem_t, elem_vn, elem_v1, elem_v2, E, nu, surface_pressure);

    elem_dofs = zeros(1, 24);
    for i = 1:4
        elem_dofs((i-1)*6 + (1:6)) = (elem_nodes(i)-1)*6 + (1:6);
    end

    K_global(elem_dofs, elem_dofs) = K_global(elem_dofs, elem_dofs) + Ke;
    F_ext(elem_dofs) = F_ext(elem_dofs) + Fe;
end

% Parse point loads
if isnumeric(loads) && numel(loads) > 1
    if size(loads, 2) == 2
        for l = 1:size(loads, 1)
            dof_idx = loads(l, 1);
            val = loads(l, 2);
            F_ext(dof_idx) = F_ext(dof_idx) + val;
        end
    elseif size(loads, 2) >= 3
        for l = 1:size(loads, 1)
            nid = loads(l, 1);
            did = loads(l, 2);
            val = loads(l, 3);
            dof_idx = (nid - 1) * 6 + did;
            F_ext(dof_idx) = F_ext(dof_idx) + val;
        end
    end
end

% Parse boundary conditions
fixed_dofs = [];
prescribed_vals = [];

if ~isempty(boundary_conditions)
    if size(boundary_conditions, 2) == 2
        fixed_dofs = boundary_conditions(:, 1);
        prescribed_vals = boundary_conditions(:, 2);
    elseif size(boundary_conditions, 2) == 3
        for b = 1:size(boundary_conditions, 1)
            nid = boundary_conditions(b, 1);
            did = boundary_conditions(b, 2);
            val = boundary_conditions(b, 3);
            dof_idx = (nid - 1) * 6 + did;
            fixed_dofs = [fixed_dofs; dof_idx];
            prescribed_vals = [prescribed_vals; val];
        end
    end
end

free_dofs = setdiff(1:total_dofs, fixed_dofs);

% Solve system
U = zeros(total_dofs, 1);
U(fixed_dofs) = prescribed_vals;

F_effective = F_ext(free_dofs) - K_global(free_dofs, fixed_dofs) * U(fixed_dofs);
U(free_dofs) = K_global(free_dofs, free_dofs) \ F_effective;

% Reaction forces
R_net = K_global * U - F_ext;

end
