function [xq, yq, zq_proj, mcmc_mse, proj_mse] = ...
    gp_experiment_run_shape(tol_thres, eps1, eps2, iter, n, ls_factor, ...
    mesh_gran, num_posteriors, desired, d, shape, fid)
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
%   xq: x-component of meshgrid for querying the surface.
%   yq: y-component of meshgrid for querying the surface.
%   zq_proj: Response value component of meshgrid for projection surface.
%   mcmc_mse: MSE of average of MCMC runs.
%   proj_mse: MSE of average of projections.

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy, x1_l, x1_h, x2_l, x2_h, x1_range, x2_range] = ...
    make_noisy_convex(n, d, shape);

%% GET SAMPLES FROM GP POSTERIOR MCMC, PROJECT EACH TO CONVEX, AND STORE.
[xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, y_nsy, ...
    ls_factor, num_posteriors);
n_gp = length(xt);

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
relative_change = (proj_mse - mcmc_mse)/mcmc_mse;
mses = [mcmc_mse proj_mse relative_change];

%% PLOT TRUE CONVEX OVER ORIGINAL DATA.
% Compute meshgrid given data. Use this in all graphs.
[xq, yq] = meshgrid(x1_l:x1_range/mesh_gran:x1_h, x2_l:x2_range/mesh_gran:x2_h);
zq_conv = griddata(xt(:, 1), xt(:, 2), ytruth_on_mcmcgrid, xq, yq);
figure; subplot(1, 3, 1);
mesh(xq, yq, zq_conv); hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title('True Convex');

%% PLOT AVG MCMC OVER ORIGINAL DATA.
zq_mcmc = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, xq, yq);
subplot(1, 3, 2);
mesh(xq, yq, zq_mcmc); hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title(sprintf('Avg MCMC (MSE = %d)', mcmc_mse));

%% PLOT AVG PROJ OVER ORIGINAL DATA.
zq_proj = griddata(xt(:, 1), xt(:, 2), avg_projs, xq, yq);
subplot(1, 3, 3);
mesh(xq, yq, zq_proj); hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title(sprintf('Avg PROJ (MSE = %d)', proj_mse));

%% SAVE FILE DATA AND FIGURE.
fprintf(fid, '%s\n', strcat(num2str(mcmc_mse, '%0.7f'), ',', ...
                            num2str(proj_mse, '%0.7f'), ',', shape));
csvwrite(sprintf('%s_x_nsy.csv', shape), x_nsy);
csvwrite(sprintf('%s_y_nsy.csv', shape), y_nsy);
csvwrite(sprintf('%s_xq.csv', shape), xq);
csvwrite(sprintf('%s_yq.csv', shape), yq);
csvwrite(sprintf('%s_zq_proj.csv', shape), zq_proj);
savefig(sprintf('%s.fig', shape));

end

