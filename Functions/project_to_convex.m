function [th_trans] = project_to_convex(n, d, x, y, eps1, eps2)
% Project x, y onto convex function.
%
% Args:
%   n: Data sample size.
%   d: Dimension of data points.
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
%   eps1: These 2 epsilons are used for convergence of the convex projection algorithm.
%   eps2: Same as above.
%
% Returns:
%   th_trans: Convex response variable.

% Standardize variables x and y.
X = zeros(n,d);
for ii=1:d
    X(:,ii) = x(:,ii) - mean(x(:,ii));  % Mean-center.
    X(:,ii) = X(:,ii)/norm(X(:,ii));    % Normalize x.
end
Y = y/norm(y);                          % Normalize y.

% Define maximum number of outer iterations.
Max_Iter = 2000;
rho0 = 1/n;	 rho = rho0;

% Find the standardized convex regression function using Sen's algorithm.
[th,xi,sq_sse, prim_feas, time_vec] = FUNCvxVanilla(X,Y,eps1,eps2, Max_Iter, rho);

% Final estimate of convex regression function at the data points.
th_trans = norm(y)*th;

end

