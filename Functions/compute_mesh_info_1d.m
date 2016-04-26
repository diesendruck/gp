function [x_l, x_h, x_range, x_grid] = compute_mesh_info_1d(x, mesh_gran)
% Given x array, mesh, boundaries, and ranges.
%
% Args:
%   x: Array of inputs.
%   mesh_gran: The number of points per dimension, to be evaluated.

%
% Returns:
%   x_l: Min of x.
%   x_h: Max of x.
%   x_range: Integer of spread of x.
%   x_grid: Mesh for x.

% Get data boundaries.
x_l = min(x); 
x_h = max(x);
x_range = x_h - x_l;

% Create surface grids.
x_grid = (x_l:x_range/mesh_gran:x_h)';

end

