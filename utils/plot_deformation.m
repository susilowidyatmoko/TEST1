function plot_deformation(nodes, elements, U, scale_factor, title_str)
% PLOT_DEFORMATION Plots initial and deformed mesh side-by-side or overlaid
%
% Inputs:
%   nodes        : N x 3 nodal coordinates
%   elements     : M x 4 element connectivity matrix
%   U            : Total displacement vector (6*N x 1)
%   scale_factor : Magnification factor for displacement
%   title_str    : Title string

if nargin < 4 || isempty(scale_factor), scale_factor = 1.0; end
if nargin < 5, title_str = 'Deformed Shell Geometry'; end

num_nodes = size(nodes, 1);
num_elem  = size(elements, 1);

% Extract translational displacements (DOFs 1, 2, 3)
disp_xyz = zeros(num_nodes, 3);
for i = 1:num_nodes
    disp_xyz(i, :) = U((i-1)*6 + (1:3))';
end

nodes_deformed = nodes + scale_factor * disp_xyz;

% Compute displacement magnitudes for color mapping
disp_mag = sqrt(sum(disp_xyz.^2, 2));

figure('Color', 'w', 'Name', title_str);
hold on; grid on; box on; axis equal;

% Plot initial mesh as wireframe
for e = 1:num_elem
    elem_nodes = elements(e, :);
    x_c = nodes(elem_nodes, 1);
    y_c = nodes(elem_nodes, 2);
    z_c = nodes(elem_nodes, 3);
    plot3([x_c; x_c(1)], [y_c; y_c(1)], [z_c; z_c(1)], ':', ...
        'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.8);
end

% Plot deformed mesh with colormap
for e = 1:num_elem
    elem_nodes = elements(e, :);
    x_d = nodes_deformed(elem_nodes, 1);
    y_d = nodes_deformed(elem_nodes, 2);
    z_d = nodes_deformed(elem_nodes, 3);
    c_d = disp_mag(elem_nodes);

    patch(x_d, y_d, z_d, c_d, 'EdgeColor', [0.1, 0.1, 0.1], 'LineWidth', 1.0);
end

colormap(jet);
cb = colorbar;
ylabel(cb, 'Displacement Magnitude');
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('%s (Scale Factor: %.1fx)', title_str, scale_factor), 'FontSize', 12, 'FontWeight', 'bold');
view(3);
hold off;

end
