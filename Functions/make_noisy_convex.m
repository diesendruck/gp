function [x_nsy, y_nsy, x1_l, x1_h, x2_l, x2_h, x1_range, x2_range] = make_noisy_convex(...
    n, d, shape)
% Make convex data, with option to normalize.
%
% Note: Each function gets a unique error variance, sig.
%
% Args:
%   n: Number of data points.
%   d: Number of dimensions of each data point.
%   shape: String indicating which convex function to use.
%
% Returns:
%   x_nsy: Matrix (nxd) of convex+noise data points.
%   y_nsy: Matrix (nx1) of noisy response values.

x_nsy = zeros(n, d);
y_nsy = zeros(n, 1);
f = zeros(n, 1);
    
if strcmp(shape, 'trough')
    x_nsy = unifrnd(-10, 10, n, d);
    sig = 1.0;
    noise = sig*randn(n, 1);

    % Compute function value.
    for ii = 1:n
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.25*x1^2;            
    end
    
elseif strcmp(shape, 'paraboloid')
    x_nsy = unifrnd(-10, 10, n, d);
    sig = 1.0;
    noise = sig*randn(n, 1);

    % Compute function value.
    for ii = 1:n
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.05*sum(x1^2 + x2^2);
    end

elseif strcmp(shape, 'wave')
    x_nsy(:, 1) = unifrnd(0, 10, n, 1);
    x_nsy(:, 2) = unifrnd(0.01, 100, n, 1);
    sig = 1.0;
    noise = sig*randn(n, 1);

    % Compute function value.
    for ii = 1:n
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.1*x1^2 - log(x2);
    end

elseif strcmp(shape, 'parabolic_cylinder');
    x_nsy = unifrnd(0, 10, n, d);
    sig = 1.0;
    noise = sig*randn(n, 1);

    % Compute function value.
    for ii = 1:n
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = -2*x1 + x2 + x1^2 - 2*x1*x2 + x2^2;
    end
    
elseif strcmp(shape, 'product');
    x_nsy = unifrnd(0, 10, n, d);
    sig = 1.0;
    noise = sig*randn(n, 1);

    % Compute function value.
    for ii = 1:n
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = 0.05*x1*x2;
    end
    
elseif strcmp(shape, 'exponential');
    x_nsy(:, 1) = unifrnd(0, 1, n, 1);
    x_nsy(:, 2) = unifrnd(0, 10, n, 1);
    sig = 1.0;
    noise = sig*randn(n, 1);

    % Compute function value.
    for ii = 1:n
        x1 = x_nsy(ii, 1); x2 = x_nsy(ii, 2);
        f(ii) = exp(x1) + x2;
    end
        
    
else
    error('Shape should be either trough or paraboloid.')
end

% Get data boundaries.
x1_l = min(x_nsy(:, 1)); x1_h = max(x_nsy(:, 1));
x2_l = min(x_nsy(:, 2)); x2_h = max(x_nsy(:, 2));
x1_range = x1_h - x1_l; x2_range = x2_h - x2_l;

% Response values = convex + noise.
y_nsy = f + noise;

end

