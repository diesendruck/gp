function [x_l, x_h, x_range, x_grid] = compute_mesh_info(x)
% Given x array, makes mesh, boundaries, and ranges.
%
% Args:
%   x: Array of inputs.
%
% Returns:
%   x1_l: Min of x.
%   x1_h: Max of x.
%   x1_range: Integer of spread of x.
%   x1_grid: Mesh for x.

% Get data boundaries.
x_l = min(x); 
x_h = max(x);
x_range = x_h - x_l;

grid_granularity = 20;
% Create surface grids.
x_grid = (x_l*1.1 : x_range/grid_granularity : x_h*1.1)';

end

