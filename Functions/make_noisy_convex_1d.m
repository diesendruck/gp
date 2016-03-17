function [x_nsy, y_nsy, x_l, x_h, x_range, x_grid, y_grid_true] = make_noisy_convex(...
    n, shape)
% Make convex data, with option to normalize.
% NOTE1: Functions appear in BOTH this file and compute_truth_from_xt_1d.m.
% NOTE2: Each function gets a unique error variance, sig.
%
% Args:
%   n: Number of data points.
%   shape: String indicating which convex function to use.
%
% Returns:
%   x_nsy: Matrix (nx1) of convex+noise data points.
%   y_nsy: Matrix (nx1) of noisy response values.
%   x_l: Min of x.
%   x_h: Max of x.
%   x_range: Max - Min, of x.
%   x_grid: Grid of x values on which to evaluate function.
%   y_grid_true: Function values at x_grid points.

d = 1;
x_nsy = zeros(n, d);
y_nsy = zeros(n, 1);
f = zeros(n, 1);

if strcmp(shape, 'parabola');
    x_nsy = unifrnd(-10, 10, n, d);
    y_nsy = zeros(n,1);
    f = zeros(n,1);
    sig = 1.0;
    noise = sig*randn(n,1);
    
    % Compute function value.
    for i = 1:n,
        f(i) = 1e-3 * (x_nsy(i,:))^4;
    end
    y_nsy = f + noise;
    
    % Compute true value over x1_grid.
    [x_l, x_h, x_range, x_grid] = compute_mesh_info_1d(x_nsy);
    y_grid_true = zeros(length(x_grid), 1);
    for i = 1:length(y_grid_true),
        y_grid_true(i) = 1e-3 * x_grid(i)^4;
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;


end

