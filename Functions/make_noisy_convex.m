function [x_nsy, y_nsy] = make_noisy_convex(n, d, shape, do_grid, ...
    data_grid_gran)
% Make convex data, with option to normalize.
% NOTE1: Functions appear in BOTH this file and compute_truth_from_xt.m.
% NOTE2: Each function gets a unique error variance, sig.
%
% Args:
%   n: Number of data points.
%   d: Number of dimensions of each data point.
%   shape: String indicating which convex function to use.
%   do_grid: Indicator for whether to generate random data, or grid data.
%   data_grid_gran: Number of points per dimension. 10 means 10x10 for d=2.
%
% Returns:
%   x_nsy: Matrix (nxd) of convex+noise data points.
%   y_nsy: Matrix (nx1) of noisy response values.

if strcmp(shape, 'trough')
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([-10 -10; 10 10], data_grid_gran);
    else
        x_nsy = unifrnd(-10, 10, n, d);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.25*x1^2;            
    end
    
elseif strcmp(shape, 'paraboloid')
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([-10 -10; 10 10], data_grid_gran);
    else
        x_nsy = unifrnd(-10, 10, n, d);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.05 * (x1^2 + x2^2);
    end

elseif strcmp(shape, 'hand')
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 1; 10 10000], data_grid_gran);
    else
        x_nsy(:, 1) = unifrnd(0, 10, n, 1);
        x_nsy(:, 2) = unifrnd(1, 10000, n, 1);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.1*x1^2 - log(x2) + 10;
    end

elseif strcmp(shape, 'parabolic_cylinder');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran);
    else
        x_nsy = unifrnd(0, 10, n, d);
    end

    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.1*(-2*x1 + x2 + x1^2 - 2*x1*x2 + x2^2);
    end
    
elseif strcmp(shape, 'wolverine');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([-10 1; 10 5], data_grid_gran);
    else
        x_nsy(:, 1) = unifrnd(-10, 10, n, 1);
        x_nsy(:, 2) = unifrnd(1, 5, n, 1);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.1 * (2*x1^2/x2 + exp(x2));
    end
    
elseif strcmp(shape, 'exponential');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 1 10], data_grid_gran);
    else
        x_nsy(:, 1) = unifrnd(0, 1, n, 1);
        x_nsy(:, 2) = unifrnd(0, 10, n, 1);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = exp(x1) + x2;
    end
        
elseif strcmp(shape, 'chair');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 -10; 10 10], data_grid_gran);
    else
        x_nsy(:, 1) = unifrnd(0, 10, n, 1);
        x_nsy(:, 2) = unifrnd(-10, 10, n, 1);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.5e-3 * (exp(x1) + x2^4);
    end

elseif strcmp(shape, 'hannah2');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([-1 -1; 1 1], data_grid_gran);
    else
        x_nsy(:, 1) = unifrnd(-1, 1, n, 1);
        x_nsy(:, 2) = unifrnd(-1, 1, n, 1);
    end
    
    l = length(x_nsy);
    y_nsy = zeros(l, 1);
    f = zeros(l, 1); 
    sig = 2.0;
    noise = sig*randn(l, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = (x1 + x2)^2;
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;

end

