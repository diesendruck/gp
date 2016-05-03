% Convex and Monotone Projections of Gaussian Processes.

% Set up constants.
verbose = 1;
do_plot = 1;
do_grid = 1;
eps1 = 1e-5;
eps2 = 1e-5;
n = 100;
d = 2;
shape = 'cm1';
data_grid_gran = 5;
ls_factor = 0.3;
num_posteriors = 50;
desired = 10;
mesh_gran = 2*data_grid_gran;
dim = sqrt(length(x_nsy));

% Make test data.
[x_nsy, ~, y_nsy, ~] = make_noisy_convex(n, d, shape, do_grid, ...
    data_grid_gran);


%% SHOW ORIGINAL CONVEX MONOTONE DATA
% Get associated data about noisy data set.
do_buffer = 0;
[~, ~, ~, ~, ~, ~, xt1, xt2, xt] = compute_mesh_info(x_nsy, mesh_gran, ...
    do_buffer);

% Plot true convex over original data.
ytruth_on_grid = compute_truth_from_xt(xt, shape);

% Also produce truth for test set.
[t1, t2] = meshgrid(unique(x_nsy(:, 1))', unique(x_nsy(:, 2))'); 
tt = [t1(:) t2(:)];
ytruth_on_test = compute_truth_from_xt(tt, shape);

if do_plot
    yq_conv = griddata(xt(:, 1), xt(:, 2), ytruth_on_grid, xt1, xt2);
    figure; subplot(2, 3, 1);
    surf(xt1, xt2, yq_conv); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title('True Convex');
end


%% COMPUTE AVERAGE FROM SAMPLES OF GP POSTERIOR MCMC, AVERAGE OF 
%% PROJECTIONS OF EACH, AND MSES OF RESPECTIVE AVERAGES.
gp_proj_start_time = tic;

% Get samples from GP posterior MCMC, project each to convex, and store.
[xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, ...
    y_nsy, ls_factor, num_posteriors, mesh_gran);
n_gp = length(xt);

n_entries = min(desired, posterior_sample_count); 
mcmcs = zeros(n_gp, n_entries);
projs = zeros(n_gp, n_entries);

for index = 1:n_entries
    if verbose
        fprintf('(%s) Projecting mcmc surface to convex, sample %d.\n', ...
            shape, index);
    end
    % Sample once from posterior, and store it as a column in mcmcs.
    y_smp = Eft_s(:, randi(posterior_sample_count));
    mcmcs(:, index) = y_smp;
    % Get convex projection of sample, and store it as a column in projs.
    y_smp_convex = project_to_convex(n_gp, d, xt, y_smp, eps1, eps2);
    projs(:, index) = y_smp_convex;
end

% Compute averages over mcmc and convex projections, respectively.
avg_mcmcs = mean(mcmcs, 2);  % Row means.
avg_projs = mean(projs, 2);

% Evaluate mcmc and proj estimates on test points.
y_mcmc_test = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, t1, t2);
y_proj_test = griddata(xt(:, 1), xt(:, 2), avg_projs, t1, t2);
mse_gp = 1/length(tt) * norm(y_mcmc_test(:) - ytruth_on_test)^2;
mse_gp_proj = 1/length(tt) * norm(y_proj_test(:) - ytruth_on_test)^2;

if do_plot
    % Plot avg MCMC over original data.
    subplot(2, 3, 2);
    yq_mcmc = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, xt1, xt2);
    surf(xt1, xt2, yq_mcmc); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('Avg GP (MSE = %d)', mse_gp));

    % Plot avg proj over original data.
    subplot(2, 3, 3); 
    yq_proj = griddata(xt(:, 1), xt(:, 2), avg_projs, xt1, xt2);
    surf(xt1, xt2, yq_proj); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('Avg GP Conv (MSE = %d)', mse_gp_proj));
end

gp_proj_time_elapsed = toc(gp_proj_start_time);


%% MONOTONE PROJECTION
f = monotone_2d(x_nsy, y_nsy);
f_mono = griddata(tt(:, 1), tt(:, 2), f(:), xt1, xt2);
mse_mono = 1/length(tt) * norm(f(:) - ytruth_on_test)^2;

subplot(2, 3, 5);
surf(xt1, xt2, f_mono); hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
title(sprintf('Monotone (MSE = %d)', mse_mono));


%% CONVEX MONOTONE PROJECTION
f = convex_monotone_2d(x_nsy, y_nsy);
f_cm = griddata(tt(:, 1), tt(:, 2), f(:), xt1, xt2);
mse_cm = 1/length(tt) * norm(f(:) - ytruth_on_test)^2;

subplot(2, 3, 6);
surf(xt1, xt2, f_cm); hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
title(sprintf('Conv+Mono (MSE = %d)', mse_cm));
