function [xt1, xt2, yq_proj, surf_mse, proj_mse, relative_change] = ...
    gp_experiment_run_shape(tol_thres, eps1, eps2, iter, n, ls_factor, ...
    mesh_gran, num_posteriors, desired, d, shape, fid, surface)
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
%   surface: String, type of surface used to fit data.
%
% Returns:
%   xt1: x1-component of meshgrid for querying the surface.
%   xt2: x2-component of meshgrid for querying the surface.
%   yq_proj: Response value component of meshgrid for projection surface.
%   surf_mse: MSE of surface (either kernel regression, MAP of GP, or 
%       average of MCMC runs.
%   proj_mse: MSE of projection (either of kernel regression, of MAP of GP,
%       or of average of projections of MCMC runs.
%   relative_change: Float of relative change of surface and projection MSEs.


%% SIMULATE RAW DATA (CONVEX + NOISE).
[x_nsy, y_nsy] = make_noisy_convex(n, d, shape);

% Get associated data about noisy data set.
[~, ~, ~, ~, ~, ~, xt1, xt2, xt] = compute_mesh_info(x_nsy, mesh_gran);

% Plot true convex over original data.
ytruth_on_grid = compute_truth_from_xt(xt, shape);
yq_conv = griddata(xt(:, 1), xt(:, 2), ytruth_on_grid, xt1, xt2);
figure; subplot(1, 3, 1);
mesh(xt1, xt2, yq_conv); hold on;
plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
title('True Convex');


%% COMPUTE KERNEL REGRESSION, ITS PROJECTION, AND RESPECTIVE MSES.
if strcmp(surface, 'kernel')
    % Select optimal bandwidth, and do kernel regression.
    h0 = rand(size(x_nsy', 1), 1)*10;
    fprintf('(%s) Optimizing for kernel regression bandwidth.', shape);
    h = Opt_Hyp_Gauss_Ker_Reg(h0, x_nsy', y_nsy');
    y_kern = zeros(size(xt1));
    for i = 1:size(xt1, 1)
        for j = 1:size(xt1, 2)
            xk = [xt1(i, j); xt2(i, j)];
            y_kern(i, j) = gaussian_kern_reg(xk, x_nsy', y_nsy', h);
        end
    end
 
    % Get convex projection of kernel regression.
    fprintf('(%s) Projecting kernel regression surface to convex.', shape);
    y_kern_proj = project_to_convex(length(xt), d, xt, y_kern(:), eps1, eps2);
    
    % Compute mses and relative change
    kern_mse = 1/length(xt) * norm(y_kern(:) - ytruth_on_grid)^2;
    proj_mse = 1/length(xt) * norm(y_kern_proj - ytruth_on_grid)^2;
    relative_change = (proj_mse - kern_mse)/kern_mse;
    
    % Save mse to return in this function.
    surf_mse = kern_mse;
    
    % Plot kernel regression and convex projection over original data.    
    subplot(1, 3, 2);
    mesh(xt1, xt2, y_kern); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Kernel Regression (MSE = %d)', kern_mse));
    
    subplot(1, 3, 3);
    yq_proj = griddata(xt(:, 1), xt(:, 2), y_kern_proj, xt1, xt2);
    mesh(xt1, xt2, yq_proj); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Projection (MSE = %d)', proj_mse));

    % Save file data and figure.
    fprintf(fid, '%s,%s,%s,%s,%s\n', surface, shape, ...
                                     num2str(kern_mse, '%0.7f'), ...
                                     num2str(proj_mse, '%0.7f'), ...
                                     num2str(relative_change, '%0.7f')); 
    savefig(sprintf('%s.fig', shape));

%% COMPUTE MAP GP SURFACE, ITS PROJECTION, AND RESPECTIVE MSES.
elseif strcmp(surface, 'map')
    disp('do map');

%% GET SAMPLES FROM GP POSTERIOR MCMC, PROJECT EACH TO CONVEX, AND STORE.
elseif strcmp(surface, 'mcmc')
    [xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, ...
        y_nsy, ls_factor, num_posteriors);
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

    %% COMPUTE AVERAGES OVER MCMC AND CONVEX PROJECTIONS, RESPECTIVELY.
    avg_mcmcs = mean(mcmcs, 2);  % Row means.
    avg_projs = mean(projs, 2);

    %% COMPUTE MSE BETWEEN TRUTH AND EACH AVERAGE.
    mcmc_mse = 1/n_gp * norm(avg_mcmcs - ytruth_on_grid)^2;
    proj_mse = 1/n_gp * norm(avg_projs - ytruth_on_grid)^2;
    relative_change = (proj_mse - mcmc_mse)/mcmc_mse;
    
    % Save mse to return in this function.
    surf_mse = mcmc_mse;

    %% PLOT AVG MCMC OVER ORIGINAL DATA.
    yq_mcmc = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, xt1, xt2);
    subplot(1, 3, 2);
    mesh(xt1, xt2, yq_mcmc); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Avg MCMC (MSE = %d)', mcmc_mse));

    %% PLOT AVG PROJ OVER ORIGINAL DATA.
    yq_proj = griddata(xt(:, 1), xt(:, 2), avg_projs, xt1, xt2);
    subplot(1, 3, 3);
    mesh(xt1, xt2, yq_proj); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 30);
    title(sprintf('Avg PROJ (MSE = %d)', proj_mse));

    %% SAVE FILE DATA AND FIGURE.
    fprintf(fid, '%s,%s,%s,%s,%s\n', surface, shape, ...
                                     num2str(mcmc_mse, '%0.7f'), ...
                                     num2str(proj_mse, '%0.7f'), ...
                                     num2str(relative_change, '%0.7f')); 
    savefig(sprintf('%s.fig', shape));
end

end

