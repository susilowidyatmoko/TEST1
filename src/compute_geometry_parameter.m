function [mu, skew_angle] = compute_geometry_parameter(elem_x)
% COMPUTE_GEOMETRY_PARAMETER Computes the adjustment parameter mu for
% geometry-dependent Gauss integration (Eq. 55 in paper).
%
% elem_x : 4 x 3 node coordinates

xi  = [ 1, -1, -1,  1];
eta = [ 1,  1, -1, -1];

xr = zeros(1, 3);
xs = zeros(1, 3);

for i = 1:4
    xr = xr + 0.25 * xi(i) * elem_x(i, :);
    xs = xs + 0.25 * eta(i) * elem_x(i, :);
end

norm_xr = norm(xr);
norm_xs = norm(xs);

if norm_xr < 1e-12 || norm_xs < 1e-12
    cos_theta = 0;
else
    % g_bar_r = xr, g_bar_s = xs evaluated at center
    % cos_theta = 1 / (|g_bar_r| |g_bar_s|) ??? No, paper Eq. (55):
    % cos_theta = (xr . xs) / (|xr| |xs|)
    cos_theta = dot(xr, xs) / (norm_xr * norm_xs);
end

mu = cos_theta^2;
skew_angle = acos(min(max(abs(cos_theta), 0), 1));

end
