% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html

%% IMPORT GPSTUFF AND SET PATHS.
if 0  % Mac version.
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/GPstuff-4.6/
    matlab_install
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/Programs/
    run_mex_commands();
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Programs')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Functions')
end

if 0  % Linux version.
    cd ~/Documents/gp/GPstuff-4.6/
    matlab_install
    cd ~/Documents/gp/Programs/
    run_mex_commands();
    cd ~/Documents/gp/
    addpath('~/Documents/gp/')
    addpath('~/Documents/gp/Programs')
    addpath('~/Documents/gp/Functions')
end

%% SET CONSTANTS.
tol_thres = 0;
eps1 = 10^-5;          % These 2 epsilons are used for convergence of the
eps2 = 10^-5;          %   convex projection algorithm.
iter = 0;              % Counter for iterations.
n = 40;                % Sample size.
ls_factor = 0.03;      % Lengthscale factor (proportion of x-range).
mesh_gran = 100;       % Number of ticks on mesh for plotting.
num_posteriors = 1020;  % Number of posterior samples to generate.
desired = 200;          % Number of posterior samples to use.
shape = 'trough';  % Shape of true convex function.
d = get_shape_dimension(shape);

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x1_l, x1_h, x2_l, x2_h, x1_range, x2_range] = make_noisy_convex(...
    n, d, shape);

%% GET SAMPLES FROM GP POSTERIOR MCMC, PROJECT EACH TO CONVEX, AND STORE.
[xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, y_nsy, ...
    ls_factor, num_posteriors);
n_gp = length(xt);

n_entries = min(desired, posterior_sample_count); 
mcmcs = zeros(n_gp, n_entries);
projs = zeros(n_gp, n_entries);

for index = 1:n_entries
    disp(sprintf('Projecting convex version of sample %d.', index));
    % Sample once from posterior, and store it as a column in mcmcs.
    y_smp = Eft_s(:, randi(posterior_sample_count));
    mcmcs(:, index) = y_smp;
    % Get convex projection of sample, and store it as a column in projs.
    y_smp_convex = project_to_convex(n_gp, d, xt, y_smp, eps1, eps2);
    projs(:, index) = y_smp_convex;
end

%% COMPUTE AVERAGES OVER MCMC AND CONVEX PROJECTIONS, RESPECTIVELY.
avg_mcmcs = mean(mcmcs, 2);  % Row means.
avg_projs = mean(projs, 2);

%% COMPUTE MSE BETWEEN TRUTH AND EACH AVERAGE.
ytruth_on_mcmcgrid = compute_truth_from_xt(xt, shape);
mcmc_mse = 1/n_gp * norm(avg_mcmcs - ytruth_on_mcmcgrid)^2;
proj_mse = 1/n_gp * norm(avg_projs - ytruth_on_mcmcgrid)^2;

%% PLOT TRUE CONVEX OVER ORIGINAL DATA.
figure;
subplot(1, 3, 1);
[xq, yq] = meshgrid(x1_l:x1_range/mesh_gran:x1_h, x2_l:x2_range/mesh_gran:x2_h);
vq = griddata(xt(:, 1), xt(:, 2), ytruth_on_mcmcgrid, xq, yq);
mesh(xq, yq, vq);
hold on
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title('True Convex');

%% PLOT AVG MCMC OVER ORIGINAL DATA.
[xq, yq] = meshgrid(x1_l:x1_range/mesh_gran:x1_h, x2_l:x2_range/mesh_gran:x2_h);
vq = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, xq, yq);
subplot(1, 3, 2);
mesh(xq, yq, vq);
hold on
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title(sprintf('Avg MCMC (MSE = %d)', mcmc_mse));

%% PLOT AVG PROJ OVER ORIGINAL DATA.
[xq, yq] = meshgrid(x1_l:x1_range/mesh_gran:x1_h, x2_l:x2_range/mesh_gran:x2_h);
vq = griddata(xt(:, 1), xt(:, 2), avg_projs, xq, yq);
subplot(1, 3, 3);
mesh(xq, yq, vq);
hold on
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title(sprintf('Avg PROJ (MSE = %d)', proj_mse));

%% SAVE FIGURE.
savefig(sprintf('%s%s.fig', [shape, datestr(clock, 0)]));


