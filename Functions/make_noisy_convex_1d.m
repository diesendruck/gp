function [x_nsy, y_nsy] = make_noisy_convex_1d(n, shape, do_grid, ...
    data_grid_gran)
% Make convex data, with option to normalize.
% NOTE1: Functions appear in BOTH this file and compute_truth_from_xt.m.
% NOTE2: Each function gets a unique error variance, sig.
%
% Args:
%   n: Number of data points.
%   shape: String indicating which convex function to use.
%   do_grid: Indicator for whether to generate random data, or grid data.
%   data_grid_gran: Number of points per dimension. 10 means 10x1 d=1.
%
% Returns:
%   x_nsy: Matrix (nx1) of convex+noise data points.
%   y_nsy: Matrix (nx1) of noisy response values.

% Setup storage for y values and function values. 
if do_grid
    len = data_grid_gran;
else
    len = n;
end
y_nsy = zeros(len, 1);
f = zeros(len, 1);

% Define x values.
if strcmp(shape, '1d-parabola');
    if do_grid
        x_nsy = linspace(-10, 10, len)';
    else
        x_nsy = unifrnd(-10, 10, len, 1);
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);
    
    % Compute function value.
    for ii = 1:len,
        f(ii) = 1e-3 * (x_nsy(ii))^4;
    end
    
elseif strcmp(shape, '1d-exponential');
    if do_grid
        x_nsy = linspace(-3, 3, len)';
    else
        x_nsy = unifrnd(-3, 3, len, 1);
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);
    
    % Compute function value.
    for ii = 1:len,
        %f(ii) = exp(x_nsy(ii)) + 0.5*x_nsy(ii)^2;
        f(ii) = exp(x_nsy(ii));
    end
    
elseif strcmp(shape, '1d-negative_entropy');
    if do_grid
        x_nsy = linspace(1, 6, len)';
    else
        x_nsy = unifrnd(1, 6, len, 1);
    end
    
    sig = 2.0;
    noise = sig*randn(len,1);
    
    % Compute function value.
    for ii = 1:len,
        f(ii) = x_nsy(ii) * log(x_nsy(ii));
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;


end

