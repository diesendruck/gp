function [x_grid, y_proj] = gp_experiment_run_shape_1d(tol_thres, eps1, ...
    eps2, iter, n, ls_factor, mesh_gran, num_posteriors, desired, d, ...
    shape, fid)
% Run gp experiment* for a particular shape.
%
% *One experiment draws many samples from the Gaussian Process posterior,
% projects each one, then compares the posterior average to the projection
% average.
%
% Args:
%   tol_thres: ?
%   eps1: These 2 epsilons are used for convergence of the convex 
%     projection algorithm.
%   eps2: Same as above.
%   iter: Counter for iterations.
%   n: Data sample size.
%   ls_factor: Lengthscale factor (proportion of x-range).
%   mesh_gran: Number of ticks on mesh for plotting.
%   num_posteriors: Number of posterior samples to generate.
%   desired: Number of posterior samples to use.
%   d: Dimension of data points.
%   shape: String, describing original convex function.
%   fid: File ID, used to store mses for average mcmc and projection.
%
% Returns:
%   x_grid: x-grid for querying the surface.
%   y_proj: Projected y value for each x in x_grid.

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x_l, x_h, x_range, x_grid, y_grid_true] = make_noisy_convex_1d(...
    n, shape);

%% GET SAMPLES FROM GP POSTERIOR MCMC, PROJECT EACH TO CONVEX, AND STORE.
[x_grid, Eft_s, posterior_sample_count] = run_gpmc_1d(x_nsy, ...
    y_nsy, ls_factor, num_posteriors);
n_gp = length(x_grid);

n_entries = min(desired, posterior_sample_count); 
mcmcs = zeros(n_gp, n_entries);
projs = zeros(n_gp, n_entries);

for index = 1:n_entries
    disp(sprintf('(%s) Projecting convex version of sample %d.', shape, ...
        index));
    % Sample once from posterior, and store it as a column in mcmcs.
    y_smp = Eft_s(:, randi(posterior_sample_count));
    mcmcs(:, index) = y_smp;
    % Get convex projection of sample, and store it as a column in projs.
    y_smp_convex = project_to_convex(n_gp, d, x_grid, y_smp, eps1, eps2);
    projs(:, index) = y_smp_convex;
end

%% COMPUTE AVERAGES OVER MCMC AND CONVEX PROJECTIONS, RESPECTIVELY.
avg_mcmcs = mean(mcmcs, 2);  % Row means.
avg_projs = mean(projs, 2);

%% COMPUTE MSE BETWEEN TRUTH AND EACH AVERAGE.
ytruth_on_mcmcgrid = compute_truth_from_xt_1d(x_grid, shape);
mcmc_mse = 1/n_gp * norm(avg_mcmcs - ytruth_on_mcmcgrid)^2;
proj_mse = 1/n_gp * norm(avg_projs - ytruth_on_mcmcgrid)^2;
relative_change = (proj_mse - mcmc_mse)/mcmc_mse;
mses = [mcmc_mse proj_mse relative_change];

%% PLOT TRUE CONVEX OVER ORIGINAL DATA.
% Compute new query grid given data. Use this in all graphs.
figure; subplot(1, 3, 1);
plot(x_grid, ytruth_on_mcmcgrid, 'k-.'); hold on;
plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
title('True Convex');

%% PLOT AVG MCMC OVER ORIGINAL DATA.
subplot(1, 3, 2);
plot(x_grid, avg_mcmcs, 'b-.'); hold on;
plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
title(sprintf('Avg MCMC (MSE = %d)', mcmc_mse));

%% PLOT AVG PROJ OVER ORIGINAL DATA.
subplot(1, 3, 3);
plot(x_grid, avg_projs, 'b-.'); hold on;
plot(x_nsy, y_nsy, 'r.', 'MarkerSize', 20);
title(sprintf('Avg PROJ (MSE = %d)', proj_mse));

%% SAVE FILE DATA AND FIGURE.
fprintf(fid, '%d,%d,%d\n', mses);
csvwrite(sprintf('Results_1d/%s_x_nsy.csv', shape), x_nsy);
csvwrite(sprintf('Results_1d/%s_y_nsy.csv', shape), y_nsy);
csvwrite(sprintf('Results_1d/%s_x_grid.csv', shape), x_grid);
csvwrite(sprintf('Results_1d/%s_avg_projs.csv', shape), avg_projs);
savefig(sprintf('Results_1d/%s.fig', shape));

end

