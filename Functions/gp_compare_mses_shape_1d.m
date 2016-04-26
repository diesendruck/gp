function [mse_gp, mse_gp_proj, mse_kern, mse_kern_proj, mse_sen, ...
    mse_cap, mse_mbcr] = gp_compare_mses_shape_1d(tol_thres, eps1, eps2,...
    iter, n, ls_factor, mesh_gran, gp_optimization, num_posteriors, ...
    desired, d, shape, fid, mbcr_burn, mbcr_tot, do_grid, data_grid_gran)
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
%   gp_optimization: "mcmc" or "map" to select GP optimization type.
%   num_posteriors: Number of posterior samples to generate.
%   desired: Number of posterior samples to use.
%   d: Dimension of data points.
%   shape: String, describing original convex function that produced data.
%   fid: File ID, used to store mses for average mcmc and projection.
%   mbcr_burn: Number of burn-in for MBCR estimate.
%   mbcr_tot: Number of total samples for MBCR estimate.
%   do_grid: Indicator for whether to generate random data, or grid data.
%   data_grid_gran: Number of points per dimension. 10 means 10x1 for d=1.
%
% Returns:
%   mse_gp: MSE of average of desired number of GP posteriors.
%   mse_gp_proj: MSE of average of projections of desired number of GP posteriors.
%   mse_kern: MSE of kernel regression over raw data.
%   mse_kern_proj: MSE of projection of kernel regression.
%   mse_sen: MSE of Sen's projection algorithm on points.
%   mse_cap: MSE of Convex Adaptive Partioning.
%   mse_mbcr: MSE of Multivariate Bayesian Covex Regression with Rvrs Jump.

% Toggle plotting on and off.
do_plot = 1;
verbose = 0;


%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy] = make_noisy_convex_1d(n, shape, do_grid, data_grid_gran);
%%
% Get associated data about noisy data set.
[~, ~, ~, x_grid] = compute_mesh_info_1d(x_nsy, mesh_gran);
n_grid = length(x_grid);

% Plot true convex over original data.
ytruth_on_grid = compute_truth_from_xt_1d(x_grid, shape);

% Also plot truth on test set.
ytruth_on_test = compute_truth_from_xt_1d(x_nsy, shape);

if do_plot
    figure; subplot(3, 3, 1);
    plot(x_grid, ytruth_on_grid, 'k-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
    title('True Convex');
end


%% COMPUTE KERNEL REGRESSION, ITS PROJECTION, AND RESPECTIVE MSES.

% Optimize kernel regression.
h0 = rand(size(x_nsy', 1), 1)*10;
h = Opt_Hyp_Gauss_Ker_Reg(h0, x_nsy', y_nsy');

% --------X GRID--------

% Kernel regression on x grid.
y_kern = zeros(size(x_grid));
for ii = 1:length(x_grid)
    y_kern(ii) = gaussian_kern_reg(x_grid(ii), x_nsy', y_nsy', h);
end

% Project to convex.
y_kern_proj = project_to_convex(length(x_grid), d, x_grid, y_kern, eps1, eps2);

% Compute mses on x grid.
%mse_kern = 1/length(x_grid) * norm(y_kern - ytruth_on_grid)^2;
%mse_kern_proj = 1/length(x_grid) * norm(y_kern_proj - ytruth_on_grid)^2;

% --------TEST GRID--------

% Kernel regression on test grid.
y_kern_test = zeros(size(x_nsy));
for ii = 1:length(x_nsy)
    y_kern_test(ii) = gaussian_kern_reg(x_nsy(ii), x_nsy', y_nsy', h);
end

% Project to convex.
y_kern_test_proj = project_to_convex(length(x_nsy), d, x_nsy, ...
    y_kern_test, eps1, eps2);

% Compute mses on test.
mse_kern = 1/length(x_nsy) * norm(y_kern_test - ytruth_on_test)^2;
mse_kern_proj = 1/length(x_nsy) * norm(y_kern_test_proj - ytruth_on_test)^2;


% Plot kernel regression and convex projection over original data.
if do_plot
    subplot(3, 3, 2);
    plot(x_grid, y_kern, 'k-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'MarkerSize', 20);
    %plot(x_nsy, y_kern_test, 'b.', 'MarkerSize', 20);
    title(sprintf('Kernel (MSE = %d)', mse_kern));
    
    subplot(3, 3, 3);
    plot(x_grid, y_kern_proj, 'k-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('Kernel Proj (MSE = %d)', mse_kern_proj));
end


%% COMPUTE GP SURFACE.
if strcmp(gp_optimization, 'mcmc');
    
    % Compute average from samples of GP Posterior MCMC, average of
    % projections of each, and MSEs of respective averages.
    
    % Get samples from GP posterior MCMC, project each to convex, and store.
    [Eft_s, posterior_sample_count] = run_gpmc_1d(x_nsy, y_nsy, ...
        ls_factor, num_posteriors, mesh_gran);

    n_entries = min(desired, posterior_sample_count); 
    mcmcs = zeros(n_grid, n_entries);
    projs = zeros(n_grid, n_entries);

    % Sample several times from posterior, and store as column in mcmcs.
    for index = 1:n_entries
        disp(sprintf('(%s) Projecting convex version of sample %d.', shape, ...
           index));
        
        % Sample posterior.
        y_smp = Eft_s(:, randi(posterior_sample_count));
        mcmcs(:, index) = y_smp;
        
        % Get convex projection of sample, and store it as a column in projs.
        y_smp_convex = project_to_convex(n_grid, d, x_grid, y_smp, eps1, eps2);
        projs(:, index) = y_smp_convex;
    end

    % Compute averages over mcmc and convex projections, respectively.
    avg_GPmcmc = mean(mcmcs, 2);  % Row means.
    avg_GPmcmc_proj = mean(projs, 2);
    
    % Redefine names for plotting (so names are consistent for either mcmc
    % vs map).
    plot_gp = avg_GPmcmc; 
    plot_proj = avg_GPmcmc_proj;

    % Compute MSE between truth and each average.
    %mse_gp = 1/n_grid * norm(avg_GPmcmc - ytruth_on_grid)^2;
    %mse_gp_proj = 1/n_grid * norm(avg_GPmcmc_projs - ytruth_on_grid)^2;
    
    % Compute MSEs for test data.
    y_mcmc_test = interp1(x_grid, avg_GPmcmc, x_nsy);
    y_proj_test = interp1(x_grid, avg_GPmcmc_proj, x_nsy);
    mse_gp = 1/length(x_nsy) * norm(y_mcmc_test - ytruth_on_test)^2;
    mse_gp_proj = 1/length(x_nsy) * norm(y_proj_test- ytruth_on_test)^2;

% elseif strcmp(gp_optimization, 'map')
%     % Get a single estimate of fitted GP.
%     [y_GPmap] = run_gp_1d(x_nsy, y_nsy, ls_factor, x_grid)
% 
%     % Project fitted GP to convex.
%     y_GPmap_proj = project_to_convex(n_grid, d, x_grid, y_GPmap, eps1, eps2);
%     
%     % Redefine names for plotting.
%     plot_gp = y_GPmap; plot_proj = y_GPmap_proj;
%     
%     % Compute MSE for map and its projection.
%     mse_gp = 1/n_grid * norm(y_GPmap - ytruth_on_grid)^2;
%     mse_gp_proj = 1/n_grid * norm(y_GPmap_proj - ytruth_on_grid)^2;
end

% Plot GP and its projection over true data.
if do_plot
    % Plot average MCMC over original data.
    subplot(3, 3, 5);
    plot(x_grid, plot_gp, 'b-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
    title(sprintf('GP (MSE = %d)', mse_gp));

    % Plot average projection over original data..
    subplot(3, 3, 6);
    plot(x_grid, plot_proj, 'b-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('GP Proj (MSE = %d)', mse_gp_proj));
end


%% COMPUTE SEN ESTIMATE AND ITS MSE.
if do_grid
    y_sen_test = project_to_convex(length(x_nsy), d, x_nsy, y_nsy, eps1, eps2);
    mse_sen = 1/length(x_nsy) * norm(y_sen_test - ytruth_on_test)^2;

    if do_plot
        % Plot avg MCMC over original data.
        subplot(3, 3, 7);
        plot(x_nsy, y_nsy, 'r.', 'MarkerSize', 20); hold on;
        plot(x_nsy, y_sen_test, 'b.', 'MarkerSize', 20);
        title(sprintf('Sen (MSE = %d)', mse_sen));
    end
else
    mse_sen = 1e10;
    if verbose
        fprintf('Skipped Sen, putting placeholder value for MSE.\n');
    end
end


%% COMPUTE CAP ESTIMATE AND ITS MSE.
[alpha, beta, K] = CAP(x_nsy, y_nsy);
y_cap = convexEval(alpha, beta, x_grid);
%mse_cap = 1/length(x_grid) * norm(y_cap - ytruth_on_grid)^2;

% Evaluate CAP estimate on test points.
y_cap_test = interp1(x_grid, y_cap, x_nsy);
mse_cap = 1/length(x_nsy) * norm(y_cap_test - ytruth_on_test)^2;

if do_plot
    % Plot average MCMC over original data.
    subplot(3, 3, 8);
    plot(x_grid, y_cap, 'b-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
    title(sprintf('CAP (MSE = %d)', mse_cap));
end


%% COMPUTE MBCR ESTIMATE AND ITS MSE.
[struct_min] = MBCR(x_nsy, y_nsy, mbcr_burn, mbcr_tot);
f_min = @(x)fMBCR(x, struct_min)
n_x_grid = length(x_grid);
y_mbcr = zeros(n_x_grid, 1);
for i = 1:n_x_grid
    y_mbcr(i) = f_min(x_grid(i));
end
%mse_mbcr = 1/length(x_grid) * norm(y_mbcr - ytruth_on_grid)^2;

% Evaluate MBCR estimate on test points.
y_mbcr_test = interp1(x_grid, y_mbcr, x_nsy);
mse_mbcr = 1/length(x_nsy) * norm(y_mbcr_test - ytruth_on_test)^2;

if do_plot
    % Plot average MCMC over original data.
    subplot(3, 3, 9);
    plot(x_grid, y_mbcr, 'b-.'); hold on;
    plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
    title(sprintf('MBCR (MSE = %d)', mse_mbcr));
end


%% ADD SUMMARY TEXT TO PLOT.
% Add text on plot to say which method did best (lowest MSE).
if do_plot
    ax = subplot(3, 3, 4);
    [~, index] = min([mse_gp mse_gp_proj mse_kern mse_kern_proj ...
                      mse_sen mse_cap mse_mbcr]);
    methods = {'mse_gp' 'mse_gp_proj' 'mse_kern' 'mse_kern_proj' ...
               'mse_sen' 'mse_cap' 'mse_mbcr'};
    min_str = strrep(char(methods(index)), '_', '\_');
    text(0, 0.5, 'Min MSE Method:', 'FontSize', 14);
    text(0, 0.3, min_str, 'FontSize', 14);
    set (ax, 'visible', 'off')

    % % Add text on plot to say which gp optimization method was used.
    % ax = subplot(3, 3, 7);
    % text(0, 0.5, 'GP Optimization:', 'FontSize', 14);
    % text(0, 0.3, gp_optimization, 'FontSize', 14);
    % set (ax, 'visible', 'off')
end


%% SAVE FILE DATA AND FIGURE.
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s\n', shape, ...
        num2str(mse_gp, '%0.7f'), ...
        num2str(mse_gp_proj, '%0.7f'), ...
        num2str(mse_kern, '%0.7f'), ...
        num2str(mse_kern_proj, '%0.7f'), ...
        num2str(mse_sen, '%0.7f'), ...
        num2str(mse_cap, '%0.7f'),...
        num2str(mse_mbcr, '%0.7f'));
                                
end
