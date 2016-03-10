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
eps1 = 10^-5;          % These 2 epsilons are used for convergence of the algo.
eps2 = 10^-5;
iter = 0;              % Counter for iterations
ls_factor = 0.03;      % Lengthscale factor (proportion of x-range)
mesh_gran = 100;       % Number of ticks on mesh for plotting.
n = 40;                % Sample size
shape = 'exponential';  % Shape of underlying convex function.
d = get_shape_dimension(shape);

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x1_l, x1_h, x2_l, x2_h, x1_range, x2_range] = make_noisy_convex(...
    n, d, shape);

%% RUN GP ON RAW DATA.
figure; subplot(2, 2, 2);
[xt_gp, y_gp] = run_gp(x_nsy, y_nsy, ls_factor, 'MAP');
title('GP MAP');
subplot(2, 2, 3);
[xt_gp, y_gp] = run_gp(x_nsy, y_nsy, ls_factor, 'MCMC');
title('GP MCMC');

%% GET CONVEX PROJECTION.
n_gp = length(xt_gp);
convex_y = project_to_convex(n_gp, d, xt_gp, y_gp, eps1, eps2);

%% PLOT TRUE CONVEX OVER ORIGINAL DATA.
subplot(2, 2, 1);
[xq, yq] = meshgrid(x1_l:x1_range/mesh_gran:x1_h, x2_l:x2_range/mesh_gran:x2_h);
ytruth_on_mcmcgrid = compute_truth_from_xt(xt_gp, shape);
vq = griddata(xt_gp(:, 1), xt_gp(:, 2), ytruth_on_mcmcgrid, xq, yq);
mesh(xq, yq, vq);
hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title('True Convex');

%% PLOT CONVEX OVER ORIGINAL GP.
[xq, yq] = meshgrid(x1_l:x1_range/mesh_gran:x1_h, x2_l:x2_range/mesh_gran:x2_h);
vq = griddata(xt_gp(:, 1), xt_gp(:, 2), convex_y, xq, yq);
subplot(2, 2, 4)
mesh(xq, yq, vq);
hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title('Convex Projection of GP MCMC');
