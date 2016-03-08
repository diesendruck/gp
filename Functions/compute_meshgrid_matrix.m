function [xt, x_range] = compute_meshgrid_matrix(x)
% Takes x matrix, and produces meshgrid matrix for surf, and boundary.
%
% Args:
%   x: Matrix of dimension nxd of inputs.
%
% Returns:
%   xt: Meshgrid matrix.
%   x_range: Biggest span across all x dimensions.

grid_granularity = 20;
x_min = min(min(x));  % Min/max across d dimensions.
x_max = max(max(x));
x_range = x_max - x_min;
x_grid = x_min*1.1 : x_range/grid_granularity : x_max*1.1;  % Create surface grid.
[xt1, xt2] = meshgrid(x_grid, x_grid);
xt = [xt1(:) xt2(:)];

end

