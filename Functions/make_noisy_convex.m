function [x_nsy, y_nsy, x_true, y_true] = make_noisy_convex(n, d, sig, shape)
% Make convex data, with option to normalize.

x_nsy = unifrnd(-10, 10, n, d);     % Design points.
y_nsy = zeros(n, 1);
f = zeros(n, 1);
noise = sig*randn(n, 1);

if d == 1
    x_true = (-10:0.2:10)';
    n_true = length(x_true);
    y_true = zeros(n_true, 1);
elseif d == 2
    [x_true1, x_true2] = meshgrid(-10:0.2:10, -10:0.2:10);
    x_true = [x_true1(:) x_true2(:)];
    n_true = length(x_true);
    y_true = zeros(n_true, 1);
end

% For each shape, store function vals over data points, and over true grid.
if strcmp(shape, 'parabola')
    for ii = 1:n
        f(ii) = x_nsy(ii)^2;
    end
    for ii = 1:n_true
        y_true(ii) = x_true(ii)^2;
    end
        
elseif strcmp(shape, 'trough')
    for ii = 1:n
        f(ii) = 0.25 * x_nsy(ii, 1)^2;            
    end
    for ii = 1:n_true
        y_true(ii) = 0.25 * x_true(ii, 1)^2;
    end
    
elseif strcmp(shape, 'paraboloid')
    for ii = 1:n,
        f(ii) = norm(x_nsy(ii,:))^2;
    end
    for ii = 1:n_true
        y_true(ii) = norm(x_true(ii, :))^2;
    end

else
    error('Shape should be either trough or paraboloid.')
end

% Response values = convex + noise.
y_nsy = f + noise;

end

