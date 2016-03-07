

s = RandStream('mt19937ar','Seed',2008)
          RandStream.setGlobalStream(s)
          

%%---- some random defintions          

tol_thres = 0;
eps1 = 10^-5;   %These 2 epsilons are used for convergence of the algo
eps2 = 10^-5;
iter = 0;   %counter for iterations


%%-----------------------

% Here is a demo example

n=100;          % Sample size
d=2;            % Dimension d
sig=0.5;        % Error variance

% Define the random variables
x = unifrnd(-1,1,n,d);      % The design points
y = zeros(n,1);
f = zeros(n,1);
nse = sig*randn(n,1);
for i = 1:n,
    f(i) = norm(x(i,:))^2;     % the function which we are performing convex regression upon
end
y = f + nse;        % Generate the response values


%%-------------------------
% standardize variables x and y
X = zeros(n,d);
for ii=1:d
    X(:,ii) = x(:,ii) - mean(x(:,ii));
    X(:,ii) = X(:,ii)/norm(X(:,ii));  % ERROR: should the normalization be wrt X, not x?
end

Y = y/norm(y);


%%-----------------------

% define the maximum number of outer iterations
Max_Iter = 2000;
rho0 = 1/n;	 rho = rho0;

% Find the standardized convex regression function 
[th,xi,sq_sse, prim_feas, time_vec] = FUNCvxVanilla(X,Y,eps1,eps2, Max_Iter, rho);

% Final estimate of convex regression function at the data points

th_trans = norm(y)*th;   

h= figure(1)
plot(cumsum(time_vec),sq_sse/n)
title('Time versus Training SSE/n')     % SSE = sum of squared errors


h= figure(2)
plot(cumsum(time_vec),prim_feas)
title('Time versus Primal Feasibility/n')





