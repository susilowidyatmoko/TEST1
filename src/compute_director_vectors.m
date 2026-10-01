function [Vn, V1, V2] = compute_director_vectors(nodes, elements)
% COMPUTE_DIRECTOR_VECTORS Computes director vectors for 4-node shell mesh
%
% Inputs:
%   nodes    : N x 3 nodal coordinates [x, y, z]
%   elements : M x 4 element connectivity matrix
%
% Outputs:
%   Vn : N x 3 normal unit director vectors
%   V1 : N x 3 first orthogonal unit vectors
%   V2 : N x 3 second orthogonal unit vectors

num_nodes = size(nodes, 1);
num_elem = size(elements, 1);

Vn = zeros(num_nodes, 3);

% Accumulate normal vectors from adjacent elements
for e = 1:num_elem
    elem_nodes = elements(e, :);
    x = nodes(elem_nodes, :);

    % Diagonal vectors for quadrilateral element
    v13 = x(3,:) - x(1,:);
    v24 = x(4,:) - x(2,:);

    % Normal vector to element mid-surface
    n_e = cross(v13, v24);
    len_n = norm(n_e);
    if len_n > 1e-12
        n_e = n_e / len_n;
    end

    for i = 1:4
        nid = elem_nodes(i);
        Vn(nid, :) = Vn(nid, :) + n_e;
    end
end

% Normalize node director vectors
for i = 1:num_nodes
    len = norm(Vn(i, :));
    if len < 1e-12
        Vn(i, :) = [0, 0, 1];
    else
        Vn(i, :) = Vn(i, :) / len;
    end
end

% Construct local orthogonal triads V1 and V2
V1 = zeros(num_nodes, 3);
V2 = zeros(num_nodes, 3);

for i = 1:num_nodes
    vn = Vn(i, :)';
    % Choose a vector not parallel to vn
    if abs(vn(1)) < 0.8 && abs(vn(2)) < 0.8
        e_ref = [1; 0; 0];
    else
        e_ref = [0; 1; 0];
    end

    v1 = cross(e_ref, vn);
    v1 = v1 / norm(v1);
    v2 = cross(vn, v1);
    v2 = v2 / norm(v2);

    V1(i, :) = v1';
    V2(i, :) = v2';
end

end
