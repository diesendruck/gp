function [th_trans] = project_to_convex(n, d, x, y, eps1, eps2, Max_Iter, rho)
% Project x, y onto convex function.

% Standardize variables x and y.
X = zeros(n,d);
for ii=1:d
    X(:,ii) = x(:,ii) - mean(x(:,ii));
    X(:,ii) = X(:,ii)/norm(X(:,ii));
end
Y = y/norm(y);
%[XT1, XT2, Z] = run_gp(X, Y, 0.1);
%scatter3(X(:,1), X(:,2), Y);

% Define maximum number of outer iterations.
Max_Iter = 2000;
rho0 = 1/n;	 rho = rho0;

% Find the standardized convex regression function.
[th,xi,sq_sse, prim_feas, time_vec] = FUNCvxVanilla(X,Y,eps1,eps2, Max_Iter, rho);

% Final estimate of convex regression function at the data points.
th_trans = norm(y)*th;

end

