% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html


s = RandStream('mt19937ar', 'Seed', 2008)
RandStream.setGlobalStream(s)       


% SET CONSTANTS -----------------------------------------------------------

tol_thres = 0;
eps1 = 10^-5;     %These 2 epsilons are used for convergence of the algo
eps2 = 10^-5;
iter = 0;         %counter for iterations
n = 20;            % Sample size
d = 2;            % Dimension d
sig = 7.0;        % Error variance

% CREATE TRAINING DATA (CONVEX + NOISE) -----------------------------------

x = unifrnd(-10, 10, n, d);
y = zeros(n, 1);
f = zeros(n, 1);
noise = sig*randn(n, 1);
for i = 1:n
    f(i) = x(i, 1)^2;  % Convex wrt first dim of x,
end                    % the other randomly varies depth of trough.
y = f + noise;
scatter3(x(:,1), x(:,2), y);

% GP REGRESSION ON TRAINING DATA ------------------------------------------
% Goal: Make GP regression data, finish with vectors x, y.

% Step 1. Train the GP.
lik = lik_gaussian('sigma2', 0.2^2);
gpcf = gpcf_sexp('lengthScale', [1.1 1.2], 'magnSigma2', 0.2^2)
pn=prior_logunif();
lik = lik_gaussian(lik, 'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gp = gp_set('lik', lik, 'cf', gpcf);

    %Should work now that X is multivariate
    %[K, C] = gp_trcov(gp, x);
    %example_x = [-1 -1 ; 0 0 ; 1 1];
    %[K, C] = gp_trcov(gp, example_x)

% Step 2. Optimize GP and get params, MAP version.
opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','iter');
gp=gp_optim(gp,x,y,'opt',opt);
[w,s] = gp_pak(gp);
disp(s), disp(exp(w))
% Step 2. Optimize GP and get params, MCMC version.  % --NOT WORKING YET--.
%[gp_rec,g,opt] = gp_mc(gp, x, y, 'nsamples', 220);
%gp_rec = thin(gp_rec, 21, 2);

% Step 3. Create surface grid.
[xt1,xt2]=meshgrid(-10.0:0.2:10.0,-10.0:0.2:10.0);
xt=[xt1(:) xt2(:)];

% Step 4. Produce surface prediction, MAP version.
[Eft_map, Varft_map] = gp_pred(gp, x, y, xt);
z = reshape(Eft_map, size(xt1));
% Step 4. Produce surface prediction, MCMC version.  % --NOT WORKING YET--.
%[Eft_s, Varft_s] = gpmc_preds(rfull, x, y, xt);
%[Eft_mc, Varft_mc] = gp_pred(gp_rec, x, y, xt);

% Step 4. Plot the surface, with the training data.
surf(xt1, xt2, z, 'FaceColor','interp', 'EdgeColor','flat', 'FaceLighting','gouraud');
axis tight; hold on;
plot3(x(:,1), x(:,2), y, 'r.', 'MarkerSize', 40);



%%-------------------------------------------------------------------------
% Start Sen, original code.

% % Define the random variables
% x = unifrnd(-1,1,n,d);      % The design points
% y = zeros(n,1);
% f = zeros(n,1);
% nse = sig*randn(n,1);
% for i = 1:n,
%     f(i) = norm(x(i,:))^2;     % the function which we are performing convex regression upon
% end
% y = f + nse;        % Generate the response values


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





