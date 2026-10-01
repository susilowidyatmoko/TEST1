function plot_convergence(mesh_sizes, errors, legend_labels, title_str)
% PLOT_CONVERGENCE Plots log-log mesh refinement convergence curves
%
% Inputs:
%   mesh_sizes    : Array or cell array of mesh densities N (e.g. [4, 8, 16, 32])
%   errors        : Matrix or cell array of relative s-norm errors Eh
%   legend_labels : Cell array of strings for legend
%   title_str     : Title string

if nargin < 4, title_str = 'Convergence Behavior'; end

figure('Color', 'w', 'Name', title_str);
hold on; grid on; box on;

markers = {'o-', 's-', '^--', 'd--', 'v-'};
colors  = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};

if iscell(errors)
    num_curves = length(errors);
    for c = 1:num_curves
        N_vals = mesh_sizes{c};
        h_vals = 1 ./ N_vals;
        err_vals = errors{c};

        m_idx = mod(c-1, length(markers)) + 1;
        c_idx = mod(c-1, length(colors)) + 1;

        loglog(h_vals, err_vals, markers{m_idx}, 'Color', colors{c_idx}, ...
            'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', colors{c_idx});
    end
else
    num_curves = size(errors, 1);
    h_vals = 1 ./ mesh_sizes;
    for c = 1:num_curves
        m_idx = mod(c-1, length(markers)) + 1;
        c_idx = mod(c-1, length(colors)) + 1;

        loglog(h_vals, errors(c, :), markers{m_idx}, 'Color', colors{c_idx}, ...
            'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', colors{c_idx});
    end
end

% Plot optimal O(h^2) reference line
h_ref = [1/32, 1/4];
err_ref = 10 * h_ref.^2;
loglog(h_ref, err_ref, 'k--', 'LineWidth', 1.5);
legend_labels{end+1} = 'Optimal O(h^2)';

xlabel('Element Size h');
ylabel('Relative Energy Error E_h');
title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
legend(legend_labels, 'Location', 'southeast', 'FontSize', 10);
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

end
