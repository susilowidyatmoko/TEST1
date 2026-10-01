function plot_mesh(nodes, elements, title_str)
% PLOT_MESH Plots original 3D finite element mesh
%
% Inputs:
%   nodes     : N x 3 nodal coordinates
%   elements  : M x 4 element connectivity matrix
%   title_str : Title string for figure

if nargin < 3, title_str = 'Finite Element Mesh'; end

figure('Color', 'w', 'Name', title_str);
hold on; grid on; box on; axis equal;

num_elem = size(elements, 1);
for e = 1:num_elem
    elem_nodes = elements(e, :);
    x_coords = nodes(elem_nodes, 1);
    y_coords = nodes(elem_nodes, 2);
    z_coords = nodes(elem_nodes, 3);

    % Draw quad patch
    patch(x_coords, y_coords, z_coords, [0.85, 0.92, 0.98], ...
        'EdgeColor', [0.1, 0.2, 0.6], 'LineWidth', 1.0, ...
        'FaceAlpha', 0.8);
end

% Plot nodes
plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.', 'MarkerSize', 8);

xlabel('X'); ylabel('Y'); zlabel('Z');
title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
view(3);
hold off;

end
