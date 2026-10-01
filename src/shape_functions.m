function [h, dh_dr, dh_ds] = shape_functions(r, s)
% SHAPE_FUNCTIONS Computes 2D bilinear shape functions for standard 4-node quad
%
% Standard isoparametric node numbering:
% Node 1: (-1, -1)
% Node 2: ( 1, -1)
% Node 3: ( 1,  1)
% Node 4: (-1,  1)

xi  = [-1,  1,  1, -1];
eta = [-1, -1,  1,  1];

h     = zeros(4, 1);
dh_dr = zeros(4, 1);
dh_ds = zeros(4, 1);

for i = 1:4
    h(i)     = 0.25 * (1 + xi(i)*r) * (1 + eta(i)*s);
    dh_dr(i) = 0.25 * xi(i)  * (1 + eta(i)*s);
    dh_ds(i) = 0.25 * eta(i) * (1 + xi(i)*r);
end

end
