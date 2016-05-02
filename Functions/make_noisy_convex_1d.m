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
f = zeros(len, 1);

% Define x values.
if strcmp(shape, '1d_parabola');
    if do_grid
        x_nsy = linspace(0, 10, len)';
    else
        x_nsy = unifrnd(0, 10, len, 1);
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);
    
    % Compute function value.
    for ii = 1:len,
        f(ii) = 0.004 * (x_nsy(ii)-7)^4;
    end
    
elseif strcmp(shape, '1d_exponential');
    if do_grid
        x_nsy = linspace(0, 10, len)';
    else
        x_nsy = unifrnd(0, 10, len, 1);
    end
    
    sig = 2.0;
    noise = sig*randn(len, 1);
    
    % Compute function value.
    for ii = 1:len,
        %f(ii) = exp(x_nsy(ii)) + 0.5*x_nsy(ii)^2;
        f(ii) = 4.5e-4 * exp(x_nsy(ii));
    end
    
elseif strcmp(shape, '1d_negative_entropy');
    if do_grid
        x_nsy = linspace(0, 10, len)';
    else
        x_nsy = unifrnd(0, 10, len, 1);
    end
    
    sig = 2.0;
    noise = sig*randn(len,1);
    
    % Compute function value.
    for ii = 1:len,
        f(ii) = 0.4 * x_nsy(ii) * log(x_nsy(ii)+1);
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;


end

