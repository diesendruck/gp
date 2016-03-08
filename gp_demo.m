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
n = 40;            % Sample size
d = 2;             % Dimension d
sig = 10.0;         % Error variance
ls_factor = 0.06;  % Lengthscale factor (proportion of x-range)

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x_true, y_true] = make_noisy_convex(n, d, sig, 'paraboloid');

% Plot true function and noisy data.
figure; subplot(2, 2, 1)
plot3(x_true(:, 1), x_true(:, 2), y_true, 'b.', 'MarkerSize', 10);
hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 40);
title('True Convex + Noisy Data');

%% RUN GP ON RAW DATA.
subplot(2, 2, 2)
[xt_gp, y_gp] = run_gp(x_nsy, y_nsy, ls_factor, 'MAP');
title('GP MAP');
subplot(2, 2, 3)
[xt_gp, y_gp] = run_gp(x_nsy, y_nsy, ls_factor, 'MCMC');
title('GP MCMC');

%% GET CONVEX PROJECTION.
n_gp = length(xt_gp);
convex_y = project_to_convex(n_gp, d, xt_gp, y_gp, eps1, eps2);

%% PLOT CONVEX OVER ORIGINAL GP.
[xq, yq] = meshgrid(-10:.2:10);
vq = griddata(xt_gp(:,1), xt_gp(:,2), convex_y, xq, yq);
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



