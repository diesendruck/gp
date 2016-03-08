function [x_nsy, y_nsy, x_true, y_true] = make_noisy_convex(n, d, sig, shape)
% Make convex data, with option to normalize.

x_nsy = unifrnd(-10, 10, n, d);     % Design points.
y_nsy = zeros(n, 1);
f = zeros(n, 1);
noise = sig*randn(n, 1);

[x_true1, x_true2] = meshgrid(-10:0.2:10, -10:0.2:10);
x_true = [x_true1(:) x_true2(:)];
y_true = zeros(size(x_true, 1), 1);

if strcmp('trough', shape)  % Convex wrt first dim of x. Other variable is random depth.
    % Evaluate noisy y response over each data point.
    for i = 1:n
        f(i) = x_nsy(i, 1)^2;            
    end
    
    % Evaluate true y response over grid.
    for j = 1:size(x_true, 1)
        y_true(j) = x_true(j, 1)^2;
    end
    
elseif strcmp('paraboloid', shape)
    % Evaluate noisy y response over each data point.
    for i = 1:n,
        f(i) = norm(x_nsy(i,:))^2;
        y_true = norm(x_true(i, :))^2;
    end
    
    % Evaluate true y response over grid.
    for j = 1:size(x_true, 1)
        y_true(j) = norm(x_true(j, :))^2;
    end
else
    error('Shape should be either trough or paraboloid.')
end

y_nsy = f + noise;  % Response values = convex + noise.

end

