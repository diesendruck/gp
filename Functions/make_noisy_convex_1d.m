function [x_nsy, y_nsy] = make_noisy_convex_1d(n, d, shape)
% Make convex data, with option to normalize.
% NOTE1: Functions appear in BOTH this file and compute_truth_from_xt.m.
% NOTE2: Each function gets a unique error variance, sig.
%
% Args:
%   n: Number of data points.
%   d: Number of dimensions of each data point (d=1 by definition).
%   shape: String indicating which convex function to use.
%
% Returns:
%   x_nsy: Matrix (nxd) of convex+noise data points.
%   y_nsy: Matrix (nx1) of noisy response values.

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
    for ii = 1:n,
        f(ii) = 1e-3 * (x_nsy(ii))^4;
    end

elseif strcmp(shape, 'exponential');
    %x_nsy = unifrnd(-3, 2, n, d);
    x_nsy = unifrnd(-3, 3, n, d);
    y_nsy = zeros(n,1);
    f = zeros(n,1);
    sig = 1.0;
    noise = sig*randn(n,1);
    
    % Compute function value.
    for ii = 1:n,
        %f(ii) = exp(x_nsy(ii)) + 0.5*x_nsy(ii)^2;
        f(ii) = exp(x_nsy(ii));
    end
    
elseif strcmp(shape, 'negative_entropy');
    x_nsy = unifrnd(0, 6, n, d);
    y_nsy = zeros(n,1);
    f = zeros(n,1);
    sig = 1.0;
    noise = sig*randn(n,1);
    
    % Compute function value.
    for ii = 1:n,
        f(ii) = x_nsy(ii) * log(x_nsy(ii));
    end
    
else
    error('Shape not recognized.')
end

% Response values = convex + noise.
y_nsy = f + noise;


end

