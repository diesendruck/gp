function [mse_gp, mse_gp_proj, mse_kern, mse_kern_proj, mse_cap, mse_mbcr] = ...
    gp_compare_mses_shape(tol_thres, eps1, eps2, iter, n, ls_factor, ...
    mesh_gran, num_posteriors, desired, d, shape, fid, mbcr_burn, mbcr_tot)
% Run gp experiment* for a particular shape.
%
% *One experiment draws many samples from the Gaussian Process posterior,
% projects each one, then compares the posterior average to the projection
% average.
%
% Args:
%   tol_thres: ?
%   eps1: These 2 epsilons are used for convergence of the convex 
%       projection algorithm.
%   eps2: Same as above.
%   iter: Counter for iterations.
%   n: Data sample size.
%   ls_factor: Lengthscale factor (proportion of x-range).
%   mesh_gran: Number of ticks on mesh for plotting.
%   num_posteriors: Number of posterior samples to generate.
%   desired: Number of posterior samples to use.
%   d: Dimension of data points.
%   shape: String, describing original convex function that produced data.
%   fid: File ID, used to store mses for average mcmc and projection.
%   mbcr_burn: Number of burn-in for MBCR estimate.
%   mbcr_tot: Number of total samples for MBCR estimate.
%
% Returns:
%   mse_gp: MSE of average of desired number of GP posteriors.
%   mse_gp_proj: MSE of average of projections of desired number of GP posteriors.
%   mse_kern: MSE of kernel regression over raw data.
%   mse_kern_proj: MSE of projection of kernel regression.
%   mse_cap: MSE of Convex Adaptive Partioning.
%   mbcr_cap: MSE of Multivariate Bayesian Covex Regression with Rvrs Jump.

% Toggle plotting on and off.
do_plot = 1;

%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy] = make_noisy_convex(n, d, shape);

% Get associated data about noisy data set.
[~, ~, ~, ~, ~, ~, xt1, xt2, xt] = compute_mesh_info(x_nsy, mesh_gran);

% Plot true convex over original data.
ytruth_on_grid = compute_truth_from_xt(xt, shape);
if do_plot
    yq_conv = griddata(xt(:, 1), xt(:, 2), ytruth_on_grid, xt1, xt2);
    figure; subplot(3, 3, 1);
    mesh(xt1, xt2, yq_conv); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title('True Convex');
end


%% COMPUTE KERNEL REGRESSION, ITS PROJECTION, AND RESPECTIVE MSES.
% Select optimal bandwidth, and do kernel regression.
h0 = rand(size(x_nsy', 1), 1)*10;
fprintf('(%s) Optimizing for kernel regression bandwidth.', shape);
h = Opt_Hyp_Gauss_Ker_Reg(h0, x_nsy', y_nsy');
y_kern = zeros(size(xt1));
for ii = 1:size(xt1, 1)
    for jj = 1:size(xt1, 2)
        xk = [xt1(ii, jj); xt2(ii, jj)];
        y_kern(ii, jj) = gaussian_kern_reg(xk, x_nsy', y_nsy', h);
    end
end

% Get convex projection of kernel regression.
fprintf('(%s) Projecting kernel regression surface to convex.', shape);
y_kern_proj = project_to_convex(length(xt), d, xt, y_kern(:), eps1, eps2);

% Compute mses and relative change
mse_kern = 1/length(xt) * norm(y_kern(:) - ytruth_on_grid)^2;
mse_kern_proj = 1/length(xt) * norm(y_kern_proj - ytruth_on_grid)^2;

% Plot kernel regression and convex projection over original data.
if do_plot
    subplot(3, 3, 2);
    mesh(xt1, xt2, y_kern); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Kernel (MSE = %d)', mse_kern));

    subplot(3, 3, 3);
    yq_proj = griddata(xt(:, 1), xt(:, 2), y_kern_proj, xt1, xt2);
    mesh(xt1, xt2, yq_proj); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Kernel Proj (MSE = %d)', mse_kern_proj));
end


%% COMPUTE AVERAGE FROM SAMPLES OF GP POSTERIOR MCMC, AVERAGE OF 
%% PROJECTIONS OF EACH, AND MSES OF RESPECTIVE AVERAGES.
% Get samples from GP posterior MCMC, project each to convex, and store.
[xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, ...
    y_nsy, ls_factor, num_posteriors, mesh_gran);
n_gp = length(xt);

n_entries = min(desired, posterior_sample_count); 
mcmcs = zeros(n_gp, n_entries);
projs = zeros(n_gp, n_entries);

for index = 1:n_entries
    fprintf('(%s) Projecting mcmc surface to convex, sample %d.\n', ...
        shape, index);
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

% Compute MSE between truth and each average.
mse_gp = 1/n_gp * norm(avg_mcmcs - ytruth_on_grid)^2;
mse_gp_proj = 1/n_gp * norm(avg_projs - ytruth_on_grid)^2;

if do_plot
    % Plot avg MCMC over original data.
    subplot(3, 3, 5);
    yq_mcmc = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, xt1, xt2);
    mesh(xt1, xt2, yq_mcmc); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Avg GP (MSE = %d)', mse_gp));

    % Plot avg proj over original data.
    subplot(3, 3, 6); 
    yq_proj = griddata(xt(:, 1), xt(:, 2), avg_projs, xt1, xt2);
    mesh(xt1, xt2, yq_proj); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Avg GP Proj (MSE = %d)', mse_gp_proj));
end


%% COMPUTE CAP ESTIMATE AND ITS MSE.
[alpha, beta, K] = CAP(x_nsy, y_nsy);
y_cap = convexEval(alpha, beta, xt);
mse_cap = 1/length(xt) * norm(y_cap - ytruth_on_grid)^2;

% Plot CAP estimate over original data.
if do_plot
    subplot(3, 3, 8); 
    yq_cap = griddata(xt(:, 1), xt(:, 2), y_cap, xt1, xt2);
    mesh(xt1, xt2, yq_cap); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('CAP (MSE = %d)', mse_cap));
end


%% COMPUTE MBCR ESTIMATE AND ITS MSE.
[struct_min] = MBCR(x_nsy, y_nsy, mbcr_burn, mbcr_tot);
f_min = @(x)fMBCR(x, struct_min)
n_xt = length(xt);
y_mbcr = zeros(n_xt, 1);
for i = 1:n_xt
    y_mbcr(i) = f_min(xt(i,:));
end
mse_mbcr = 1/length(xt) * norm(y_mbcr - ytruth_on_grid)^2;

% Plot MBCR estimate over original data.
if do_plot
    subplot(3, 3, 9); 
    yq_mbcr = griddata(xt(:, 1), xt(:, 2), y_mbcr, xt1, xt2);
    mesh(xt1, xt2, yq_mbcr); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('MBCR (MSE = %d)', mse_mbcr));
end


%% ADD SUMMARY TEXT TO PLOT.
% Add text on plot to say which method did best (lowest MSE).
ax = subplot(3, 3, 4);
[~, index] = min([mse_gp mse_gp_proj mse_kern mse_kern_proj mse_cap mse_mbcr]);
methods = {'mse_gp' 'mse_gp_proj' 'mse_kern' 'mse_kern_proj' 'mse_cap' 'mse_mbcr'};
min_str = strrep(char(methods(index)), '_', '\_');
text(0, 0.5, 'Min MSE Method:', 'FontSize', 14);
text(0, 0.3, min_str, 'FontSize', 14);
set (ax, 'visible', 'off')

% % Add text on plot to say which gp optimization method was used.
% ax = subplot(3, 3, 7);
% text(0, 0.5, 'GP Optimization:', 'FontSize', 14);
% text(0, 0.3, gp_optimization, 'FontSize', 14);
% set (ax, 'visible', 'off')


%% SAVE FILE DATA AND FIGURE.
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n', shape, num2str(mse_gp, '%0.7f'), ...
                                    num2str(mse_gp_proj, '%0.7f'), ...
                                    num2str(mse_kern, '%0.7f'), ...
                                    num2str(mse_kern_proj, '%0.7f'), ...
                                    num2str(mse_cap, '%0.7f'),...
                                    num2str(mse_mbcr, '%0.7f')); 

end

