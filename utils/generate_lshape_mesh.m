function [nodes, elements] = generate_lshape_mesh(L, N)
% GENERATE_LSHAPE_MESH Generates mesh for L-shaped structure (Fig. 10)
% Composed of 3 square blocks of size L/2 x L/2:
% Block 1: [0, L/2] x [0, L/2]
% Block 2: [L/2, L] x [0, L/2]
% Block 3: [0, L/2] x [L/2, L]

N2 = N / 2;
x_pts = linspace(0, L, N + 1);
y_pts = linspace(0, L, N + 1);

nodes = [];
node_map = zeros(N+1, N+1);

% Create nodes (only in L-shape domain)
node_idx = 1;
for j = 1:N+1
    for i = 1:N+1
        x = x_pts(i);
        y = y_pts(j);

        % Skip top-right quadrant (x > L/2 and y > L/2)
        if i > N2 + 1 && j > N2 + 1
            continue;
        end

        nodes = [nodes; x, y, 0];
        node_map(i, j) = node_idx;
        node_idx = node_idx + 1;
    end
end

% Create element connectivity
elements = [];
for j = 1:N
    for i = 1:N
        if i > N2 && j > N2
            continue; % skip empty top-right quadrant
        end

        n1 = node_map(i, j);
        n2 = node_map(i+1, j);
        n3 = node_map(i+1, j+1);
        n4 = node_map(i, j+1);

        elements = [elements; n1, n2, n3, n4];
    end
end

end
