function [x1_l, x1_h, x2_l, x2_h, x1_range, x2_range, xt1, xt2, xt] = ...
    compute_mesh_info(x, mesh_gran)
% Given x matrix, makes meshgrid matrix, boundaries, and ranges.
%
% Args:
%   x: Matrix of dimension nxd of inputs.
%   mesh_gran: The number of points per dimension, to be evaluated.
%
% Returns:
%   x1_l: Min of first dimension of x.
%   x1_h: Max of first dimension of x.
%   x2_l: Min of second dimension of x.
%   x2_h: Max of second dimension of x.
%   x1_range: Integer of spread of first dimension of x.
%   x2_range: Integer of spread of second dimension of x.
%   xt1: First dimension component for meshgrid matrix.
%   xt2: Second dimension component for meshgrid matrix.
%   xt: Full meshgrid matrix.

% Get data boundaries.
x1_l = min(x(:, 1)); x1_h = max(x(:, 1));
x2_l = min(x(:, 2)); x2_h = max(x(:, 2));
x1_range = x1_h - x1_l; x2_range = x2_h - x2_l;

x1_grid = linspace(x1_l, x1_h, mesh_gran);  % Create surface grids.
x2_grid = linspace(x2_l, x2_h, mesh_gran);  % Create surface grids.
[xt1, xt2] = meshgrid(x1_grid, x2_grid);
xt = [xt1(:) xt2(:)];

end

