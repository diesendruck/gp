function [x, y] = make_noisy_convex(n, d, sig, shape)
% Make convex data, with option to normalize.

x = unifrnd(-10, 10, n, d);     % Design points.
y = zeros(n, 1);
f = zeros(n, 1);
noise = sig*randn(n, 1);

if strcmp('trough', shape)    
    for i = 1:n
        f(i) = x(i, 1)^2;           % Convex wrt first dim of x,
    end                             % the other randomly varies depth of trough.
elseif strcmp('paraboloid', shape)
    for i = 1:n,
        f(i) = norm(x(i,:))^2;      % Convex function (upper half of paraboloid).
    end
else
    error('Shape should be either trough or paraboloid.')
end
y = f + noise;                      % Response values = convex + noise.

end

