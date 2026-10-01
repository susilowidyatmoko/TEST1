function [nodes, elements] = generate_quad_mesh(x_span, y_span, Nx, Ny, distortion_type)
% GENERATE_QUAD_MESH Generates regular or distorted 4-node quad mesh on rectangular domain
%
% Standard isoparametric node numbering (counter-clockwise starting at (-1, -1)):
% 1: (x_i, y_j)
% 2: (x_i+1, y_j)
% 3: (x_i+1, y_j+1)
% 4: (x_i, y_j+1)

if nargin < 5, distortion_type = 'none'; end

x_min = x_span(1); x_max = x_span(2);
y_min = y_span(1); y_max = y_span(2);

% Grid coordinates
if strcmp(distortion_type, 'none')
    u = linspace(0, 1, Nx + 1);
    v = linspace(0, 1, Ny + 1);
else
    % Distorted ratio 1 : 2 : 3 : ... : N (paper pattern)
    ratios_x = 1:Nx;
    u = [0, cumsum(ratios_x) / sum(ratios_x)];
    ratios_y = 1:Ny;
    v = [0, cumsum(ratios_y) / sum(ratios_y)];
end

[U_grid, V_grid] = meshgrid(u, v);

% Generate nodes
nodes = zeros((Nx+1)*(Ny+1), 3);
idx = 1;
for j = 1:Ny+1
    for i = 1:Nx+1
        u_val = U_grid(j, i);
        v_val = V_grid(j, i);

        if strcmp(distortion_type, 'pattern2')
            if u_val > 0 && u_val < 1 && v_val > 0 && v_val < 1
                u_val = u_val + 0.05 * sin(pi * u_val) * sin(pi * v_val);
                v_val = v_val + 0.05 * sin(pi * u_val) * sin(pi * v_val);
            end
        end

        nodes(idx, :) = [x_min + u_val*(x_max - x_min), y_min + v_val*(y_max - y_min), 0];
        idx = idx + 1;
    end
end

% Generate element connectivity
elements = zeros(Nx*Ny, 4);
elem_idx = 1;
for j = 1:Ny
    for i = 1:Nx
        n1 = (j-1)*(Nx+1) + i;       % (-1, -1)
        n2 = (j-1)*(Nx+1) + i + 1;   % ( 1, -1)
        n3 = j*(Nx+1)     + i + 1;   % ( 1,  1)
        n4 = j*(Nx+1)     + i;       % (-1,  1)

        elements(elem_idx, :) = [n1, n2, n3, n4];
        elem_idx = elem_idx + 1;
    end
end

end
