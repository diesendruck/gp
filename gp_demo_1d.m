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
eps1 = 10^-5;          % These epsilons are used for algo convergence.
eps2 = 10^-5;
iter = 0;              % Counter for iterations.
ls_factor = 0.01;      % Lengthscale factor (proportion of x-range).
n = 20;                % Sample size.
d = 1;                 % Dimension of data.
shape = 'parabola';    % True convex function.

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x_l, x_h, x_range, x_grid, y_grid_true] = make_noisy_convex_1d(n, shape);

% Plot true function and noisy data.
figure; subplot(2, 2, 1);
plot(x_grid, y_grid_true, 'k-.'); hold on;
plot(x_nsy, y_nsy, 'r.', 'markers', 20);
xlim([min(x_grid) max(x_grid)]);
ylim([min(y_nsy)*1.1 max(y_nsy)*1.1]);
title('True Convex + Noisy Data');

%% RUN GP ON RAW DATA.
subplot(2, 2, 2);
[y_gp_map] = run_gp_1d(x_nsy, y_nsy, ls_factor, 'MAP', x_grid);
title('GP MAP');
subplot(2, 2, 3)
[y_gp_mcmc] = run_gp_1d(x_nsy, y_nsy, ls_factor, 'MCMC', x_grid);
title('GP MCMC');

%% GET CONVEX PROJECTION OF MCMC.
n = length(x_grid);
convex_y = project_to_convex(n, d, x_grid, y_gp_mcmc, eps1, eps2);

%% PLOT CONVEX OVER ORIGINAL GP.
subplot(2, 2, 4)
plot(x_grid, convex_y, 'b-.', 'Markers', 10);
hold on
plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
xlim([min(x_grid) max(x_grid)]);
title('Convex Projection of GP MCMC');
