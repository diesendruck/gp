function [x_nsy, x_nsy_jit, y_nsy, y_nsy_jit] = make_noisy_convex(n, d, ...
    shape, do_grid, ...
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
%   x_nsy_jit: Same as x_nsy, but jittered.
%   y_nsy: Matrix (nx1) of noisy response values.
%   y_nsy_jit: Same as y_nsy, but for jittered inputs.

% Setup storage for y values and function values. Given data_grid_gran, can
% find number of output points: e.g. grid gran of 5 means 5x5 or 25 points.
len = data_grid_gran^2;
y_nsy = zeros(len, 1);
y_nsy_jit = zeros(len, 1);
f = zeros(len, 1);
f_jit = zeros(len, 1);
do_buffer = 0;

if strcmp(shape, 'trough')
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([-10 -10; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];
    else
        x_nsy = unifrnd(-10, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1);
        x1_jit = x_nsy_jit(ii, 1);
        f(ii) = 0.25*x1^2;
        f_jit(ii) = 0.25*x1_jit^2;
    end
    
elseif strcmp(shape, 'paraboloid')
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([-10 -10; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];
    else
        x_nsy = unifrnd(-10, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.05 * (x1^2 + x2^2);
        f_jit(ii) = 0.05 * (x1_jit^2 + x2_jit^2);
    end

elseif strcmp(shape, 'hand')
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([0 1; 10 10000], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];
    else
        x_nsy(:, 1) = unifrnd(0, 10, n, 1);
        x_nsy(:, 2) = unifrnd(1, 10000, n, 1);
        x_nsy_jit = x_nsy;
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 1e-5*x1^6 - log(x2) + 10;
        f_jit(ii) = 1e-5*x1_jit^6 - log(x2_jit) + 10;
    end

elseif strcmp(shape, 'parabolic_cylinder');
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];        
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end

    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.1*(-2*x1 + x2 + x1^2 - 2*x1*x2 + x2^2);
        f_jit(ii) = 0.1*(-2*x1_jit + x2_jit + x1_jit^2 - ...
                         2*x1_jit*x2_jit + x2_jit^2);
    end
    
elseif strcmp(shape, 'wolverine');
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([-10 1; 10 5], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];
    else
        x_nsy(:, 1) = unifrnd(-10, 10, n, 1);
        x_nsy(:, 2) = unifrnd(1, 5, n, 1);
        x_nsy_jit = x_nsy;
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.1 * (2*x1^2/x2 + exp(x2));
        f_jit(ii) = 0.1 * (2*x1_jit^2/x2_jit + exp(x2_jit));
    end
    
elseif strcmp(shape, 'exponential');
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 1 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];
    else
        x_nsy(:, 1) = unifrnd(0, 1, n, 1);
        x_nsy(:, 2) = unifrnd(0, 10, n, 1);
        x_nsy_jit = x_nsy;
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = exp(x1) + x2;
        f_jit(ii) = exp(x1_jit) + x2_jit;
    end
        
elseif strcmp(shape, 'chair');
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([0 -10; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];    
    else
        x_nsy(:, 1) = unifrnd(0, 10, n, 1);
        x_nsy(:, 2) = unifrnd(-10, 10, n, 1);
        x_nsy_jit = x_nsy;
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.5e-3 * (exp(x1) + x2^4);
        f_jit(ii) = 0.5e-3 * (exp(x1_jit) + x2_jit^4);
    end

elseif strcmp(shape, 'hannah2');
    if do_grid
        [~, ~, ~, ~, x1_range, x2_range, ~, ~, x_nsy] = ...
            compute_mesh_info([-1 -1; 1 1], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        x_nsy_jit = x_nsy + [(1e-6)*x1_range*randn(length(x_nsy), 1) ...
                             (1e-6)*x2_range*randn(length(x_nsy), 1)];
    else
        x_nsy(:, 1) = unifrnd(-1, 1, n, 1);
        x_nsy(:, 2) = unifrnd(-1, 1, n, 1);
        x_nsy_jit = x_nsy;
    end

    sig = 0.5;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:length(x_nsy)
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = (x1 + x2)^2;
        f_jit(ii) = (x1_jit + x2_jit)^2;
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;
y_nsy_jit = f_jit + noise;

end

