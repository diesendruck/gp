% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html

%% IMPORT GPSTUFF AND SET PATHS.
if 0
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/GPstuff-4.6/
    matlab_install
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Programs')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Functions')
end

%% SET CONSTANTS.
tol_thres = 0;
eps1 = 10^-5;      % These 2 epsilons are used for convergence of the algo.
eps2 = 10^-5;
iter = 0;          % Counter for iterations
n = 20;            % Sample size
d = 2;             % Dimension d
sig = 5.0;         % Error variance
ls_factor = 0.06;  % Lengthscale factor (proportion of x-range)

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x_true, y_true] = make_noisy_convex(n, d, sig, 'paraboloid');

% Plot true function and noisy data.
figure; subplot(1, 3, 1)
plot3(x_true(:, 1), x_true(:, 2), y_true, 'b.', 'MarkerSize', 10);
hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 40);
title('True Convex + Noisy Data');

%% GET SAMPLES FROM GP POSTERIOR MCMC, PROJECT EACH TO CONVEX, AND STORE.
[xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, y_nsy, ...
    ls_factor);
n_gp = length(xt);
gp_dim = [length(xt1) length(xt1)];

% Declare how many posterior samples to use.
desired = 10;
n_entries = min(desired, posterior_sample_count); 
mcmcs = zeros(n_gp, n_entries);
projs = zeros(n_gp, n_entries);

for index = 1:n_entries
    % Sample one from posterior, and store it.
    y_smp = Eft_s(:, randi(posterior_sample_count));
    mcmcs(:, index) = y_smp;
    % Get convex projection of sample.
    y_smp_convex = project_to_convex(n_gp, d, xt, y_smp, eps1, eps2);
    projs(:, index) = y_smp_convex;
end




%% COMPUTE AVERAGES OVER MCMC AND CONVEX PROJECTIONS, RESPECTIVELY.


%% COMPUTE MSE BETWEEN TRUTH AND EACH AVERAGE.


%% PLOT CONVEX OVER ORIGINAL GP.
[xq, yq] = meshgrid(-10:.2:10);
vq = griddata(xt(:,1), xt(:,2), y_smp_convex, xq, yq);
subplot(2, 2, 4)
mesh(xq,yq,vq);
hold on
plot3(x_nsy(:,1), x_nsy(:,2), y_nsy, 'r.', 'MarkerSize', 40);
title('Convex Projection of GP MCMC');


%% EXTRAS.
% DIAGNOSTICS FROM SEN, on outputs of project_to_convex that aren't
% currently being returned.
% h= figure(1)
% plot(cumsum(time_vec),sq_sse/n)
% title('Time versus Training SSE/n')     % SSE = sum of squared errors
% 
% h= figure(2)
% plot(cumsum(time_vec),prim_feas)
% title('Time versus Primal Feasibility/n')



