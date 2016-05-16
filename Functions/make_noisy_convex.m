function [x_nsy, x_nsy_jit, y_nsy, y_nsy_jit] = make_noisy_convex(n, d, ...
    shape, do_grid, data_grid_gran)
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
%   x_nsy: Matrix of convex+noise data points.
%   x_nsy_jit: Same as x_nsy, but jittered.
%   y_nsy: Matrix (nx1) of noisy response values.
%   y_nsy_jit: Same as y_nsy, but for jittered inputs.

% Toggle plotting. Set to 0 when running full test in gp_compare_mses.m.
do_plot = 0;

% Setup storage for y values and function values. Given data_grid_gran, can
% find number of output points: e.g. grid gran of 5 means 5x5 or 25 points.
if do_grid
    len = data_grid_gran^2;
else
    len = n;
end
y_nsy = zeros(len, 1);
y_nsy_jit = zeros(len, 1);
f = zeros(len, 1);
f_jit = zeros(len, 1);
do_buffer = 0;

if strcmp(shape, 'trough')
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1);
        x1_jit = x_nsy_jit(ii, 1);
        f(ii) = 0.2*(x1-3)^2;
        f_jit(ii) = 0.2*(x1_jit-3)^2;
    end
    
elseif strcmp(shape, 'paraboloid')
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.2 * ((x1-5)^2 + (x2-5)^2);
        f_jit(ii) = 0.2 * ((x1_jit-5)^2 + (x2_jit-5)^2);
    end

elseif strcmp(shape, 'hand')
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 1e-5*x1^6 - log(100*x2+1000) + 8;
        f_jit(ii) = 1e-5*x1_jit^6 - log(100*x2_jit+1000) + 8;
    end

elseif strcmp(shape, 'parabolic_cylinder');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);        
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end

    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.1*(-2*x1 + x2 + x1^2 - 2*x1*x2 + x2^2);
        f_jit(ii) = 0.1*(-2*x1_jit + x2_jit + x1_jit^2 - ...
                         2*x1_jit*x2_jit + x2_jit^2);
    end
    
elseif strcmp(shape, 'wolverine');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.1 * (2*(x1-5)^2/(x2+1) + exp(2*(x2+1)/5));
        f_jit(ii) = 0.1 * (2*(x1_jit-5)^2/(x2_jit+1) + exp(2*(x2_jit+1)/5));
    end
    
elseif strcmp(shape, 'exponential');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.3*(exp(x1/3) + 0.3*x2);
        f_jit(ii) = 0.3*(exp(x1_jit/3) + 0.3*x2_jit);
    end
        
elseif strcmp(shape, 'chair');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);    
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end
    
    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 1e-3 * (exp(x1/1.1) + 2*(x2-5)^4);
        f_jit(ii) = 1e-3 * (exp(x1_jit/1.1) + 2*(x2_jit-5)^4);
    end

elseif strcmp(shape, 'hannah2');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([-1 -1; 1 1], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy(:, 1) = unifrnd(-1, 1, n, 1);
        x_nsy(:, 2) = unifrnd(-1, 1, n, 1);
        x_nsy_jit = x_nsy;
    end

    sig = 0.5;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = (x1 + x2)^2;
        f_jit(ii) = (x1_jit + x2_jit)^2;
    end
    
elseif strcmp(shape, 'cm1');
    if do_grid
        [~, ~, ~, ~, ~, ~, ~, ~, x_nsy] = ...
            compute_mesh_info([0 0; 10 10], data_grid_gran, do_buffer);
        % Jitter the data. For CAP method.
        [x_nsy_jit] = jitter(x_nsy);
    else
        x_nsy = unifrnd(0, 10, n, d);
        x_nsy_jit = x_nsy;
    end

    sig = 3.0;
    noise = sig*randn(len, 1);

    % Compute function value.
    for ii = 1:len
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        x1_jit = x_nsy_jit(ii, 1); x2_jit = x_nsy_jit(ii, 2);
        f(ii) = 0.025 * (x1 + x2)^2;
        f_jit(ii) = 0.025 * (x1_jit + x2_jit)^2;
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;
y_nsy_jit = f_jit + noise;

if do_plot
    %p3(x_nsy, y_nsy);
    % Plot the true surface, with noisy points on top.
    do_buffer = 0;
    [~, ~, ~, ~, ~, ~, xt1, xt2, xt] = compute_mesh_info(x_nsy, mesh_gran, ...
        do_buffer);
    yq = griddata(x_nsy(:, 1), x_nsy(:, 2), f, xt1, xt2);
    mesh(xt1, xt2, yq); hold on;
    plot3(x_nsy_jit(:, 1), x_nsy_jit(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(strrep(shape, '_', '\_'));
end

end
