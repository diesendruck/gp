% Test kernel, kernel projection, and Sen methods on sample data sets.

% Constants.
start_time = tic;
n_data = 100;  % Number of data points.
d = 2;
do_plot = 1;
do_grid = 1;
data_grid_gran = 10;
mesh_gran = data_grid_gran*2;
num_posteriors = 2000;   % For run_gpmc, number of posteriors to fetch.
desired = 50;           % For run_gpmc, number of posteriors to use.
ls_factor = 0.5;
eps1 = 1e-4;              
eps2 = 1e-4;
verbose = 1;
platform = 'mac';       % Choose appropriate platform.

% Choose shape. This script runs for one shape.
shape = 'exponential';

% Import data.                           
filename = sprintf('data/50data_samples_%s.csv', shape);
xydat = importdata(filename);
x_nsy = xydat(:, [1:2]);
y_samples = xydat(:, [3:size(xydat, 2)]);
num_samples = size(y_samples, 2);

% Setup directories and email.
setup_directories_and_email(platform);

% Set up storage variable for each method, reporting RMSEs.
store_rmses_kern = zeros(num_samples, 1);
store_rmses_kern_mono = zeros(num_samples, 1);
store_rmses_gp = zeros(num_samples, 1);
store_rmses_gp_mono = zeros(num_samples, 1);


%% HELPER FUNCTIONS FROM GP_COMPARE_MSES_SHAPE, TO SET UP MESHES.
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


%% DO GLOBAL RUNS
% For each set of y outputs, produce the RMSEs for kernel, monotone
% projection of kernel, GP, and monotone projection of GP.

for i = 1:num_samples
    % Get y sample for this iteration.
    y_nsy = y_samples(:, i);
    
    %% KERNEL REGRESSION.
    h0 = rand(size(x_nsy', 1), 1)*10;
    if verbose
        fprintf('(%s) Optimizing for kernel regression bandwidth.\n', shape);
    end
    h = Opt_Hyp_Gauss_Ker_Reg(h0, x_nsy', y_nsy');
    y_kern = zeros(size(xt1));
    for ii = 1:size(xt1, 1)
        for jj = 1:size(xt1, 2)
            xk = [xt1(ii, jj); xt2(ii, jj)];
            y_kern(ii, jj) = gaussian_kern_reg(xk, x_nsy', y_nsy', h);
        end
    end
    y_kern = y_kern(:);

    % Project to monotone.
    [y_kern_mono, kern_m_iter] = monotone_2d(xt, y_kern, shape);
    % Correct the shape after kernel regression.
    y_kern_mono = y_kern_mono(:);

    % Evaluate kernel on test points.
    y_kern_test = griddata(xt(:, 1), xt(:, 2), y_kern, t1, t2);
    y_kern_mono_test = griddata(xt(:, 1), xt(:, 2), y_kern_mono, t1, t2);

    % Compute RMSE for kernel and kernel monotone.
    rmse_kern = sqrt(1/length(tt) * norm(y_kern_test(:) - ytruth_on_test)^2);
    rmse_kern_mono = sqrt(1/length(tt) * norm(y_kern_mono_test(:) - ytruth_on_test)^2);

    % Store results.
    store_rmses_kern(i, :) = rmse_kern;
    store_rmses_kern_mono(i, :) = rmse_kern_mono;
    
    if do_plot
        % Plot kernel data.
        figure
        subplot(2, 3, 1)
        mesh(xt1, xt2, reshape(ytruth_on_grid, mesh_gran, mesh_gran)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zl = zlim;
        title(sprintf('Monotone Truth'))

        subplot(2, 3, 2)
        mesh(xt1, xt2, reshape(y_kern, mesh_gran, mesh_gran)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('Kernel Unconstrained (RMSE = %s)', num2str(rmse_kern, '%0.3f')));

        subplot(2, 3, 3)
        mesh(xt1, xt2, reshape(y_kern_mono, mesh_gran, mesh_gran)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('Kernel Monotone (RMSE = %s)', num2str(rmse_kern_mono, '%0.3f')));
    end
    
    
    %% GP REGRESSION.
    % Get samples from GP posterior MCMC, project each to convex, and store.
    [xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, ...
        y_nsy, ls_factor, num_posteriors, mesh_gran);
    n_gp = length(xt);

    n_entries = min(desired, posterior_sample_count); 
    mcmcs = zeros(n_gp, n_entries);

    for index = 1:n_entries
        if verbose
            fprintf('(%s) Projecting mcmc surface to convex, sample %d.\n', ...
                shape, index);
        end
        % Sample once from posterior, and store it as a column in mcmcs.
        y_smp = Eft_s(:, randi(posterior_sample_count));
        mcmcs(:, index) = y_smp;
        
    end

    % Compute averages over mcmc and convex projections, respectively.
    avg_mcmcs = mean(mcmcs, 2);  % Row means.
    y_gp = avg_mcmcs;

    % Project to monotone.
    [y_gp_mono, gp_m_iter] = monotone_2d(xt, y_gp, shape);
    
    % Evaluate gp on test points.
    y_gp_test = griddata(xt(:, 1), xt(:, 2), y_gp, t1, t2);
    y_gp_mono_test = griddata(xt(:, 1), xt(:, 2), y_gp_mono(:), t1, t2);
    
    % Compute RMSE for GP and GP monotone.
    rmse_gp = sqrt(1/length(tt) * norm(y_gp_test(:) - ytruth_on_test)^2);
    rmse_gp_mono = sqrt(1/length(tt) * norm(y_gp_mono_test(:) - ytruth_on_test)^2);

    % Store results.
    store_rmses_gp(i, :) = rmse_gp;
    store_rmses_gp_mono(i, :) = rmse_gp_mono;
    
    if do_plot
        % Plot avgGP data.
        subplot(2, 3, 5)
        mesh(xt1, xt2, reshape(y_gp, mesh_gran, mesh_gran)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('GP Unconstrained (RMSE = %s)', num2str(rmse_gp, '%0.3f')));

        subplot(2, 3, 6)
        mesh(xt1, xt2, reshape(y_gp_mono, mesh_gran, mesh_gran)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('GP Monotone (RMSE = %s)', num2str(rmse_gp_mono, '%0.3f')));
    end
    
end

toc(start_time);

% Save results to file.
results = [store_rmses_gp store_rmses_gp_mono store_rmses_kern ...
           store_rmses_kern_mono];
csvwrite(sprintf('data/rmses_gp_gp_m_kern_kern_m_50data_%s.csv', shape), results);

if strcmp(platform, 'mac')
    sendmail('momod@utexas.edu', ...
        sprintf('RESULTS mac: shively_%s_rmses', shape), ...
        sprintf('Total time: %s', num2str(toc(start_time), '%0.2f')), ...
        {strcat('/Users/mauricediesendruck/Google Drive/', ...
                '0-LIZHEN RESEARCH/gp/data/', ...
                sprintf('rmses_gp_gp_m_kern_kern_m_50data_%s.csv', shape))});       
elseif strcmp(platform, 'linux')
    sendmail('momod@utexas.edu', ...
        sprintf('RESULTS linux: shively_%s_rmses', shape), ...
        sprintf('Total time: %s', num2str(toc(start_time), '%0.2f')), ...
        {strcat('/home/momod/Documents/gp/data/', ...
                sprintf('rmses_gp_gp_m_kern_kern_m_50data_%s.csv', shape))});
end

    
    