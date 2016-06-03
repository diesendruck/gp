% Test kernel, kernel projection, and Sen methods on sample data sets.

% Constants.
start_time = tic;
n_data = 100;  % Number of data points.
d = 2;
do_plot = 1;
do_grid = 1;
data_grid_gran = 10;
dim = data_grid_gran;
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
filename = 'data/UnconstrainedFunctionEstimates_50Datasets_Neq100_Exponential.dat';
xydat = importdata(filename);
x_nsy = xydat(:, [1:2]);
y_samples = xydat(:, [3:size(xydat, 2)]);
num_samples = size(y_samples, 2);

% Setup directories and email.
setup_directories_and_email(platform);

% Set up storage variable for each method, reporting RMSEs.
store_rmses = zeros(num_samples, 1);
store_rmses_mono = zeros(num_samples, 1);
store_rmses_conv = zeros(num_samples, 1);
store_rmses_cm = zeros(num_samples, 1);


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
    %% Get y sample for this iteration.
    y_nsy = y_samples(:, i);

    % Compute RMSE for spline, and store results.
    rmse = sqrt(1/length(tt) * norm(y_nsy - ytruth_on_test)^2);
    store_rmses(i, :) = rmse;
    
    %% Project to monotone.
    [y_mono, m_iter] = monotone_2d(tt, y_nsy, shape);
    y_mono = y_mono(:);

    % Compute RMSE for spline monotone, and store results.
    rmse_mono = sqrt(1/length(tt) * norm(y_mono - ytruth_on_test)^2);
    store_rmses_mono(i, :) = rmse_mono;
    
    
    %% Project spline to convex.
    y_conv = project_to_convex(length(tt), d, tt, y_nsy, eps1, eps2);

    % Compute RMSE for spline convex, and store results.
    rmse_conv = sqrt(1/length(tt) * norm(y_conv - ytruth_on_test)^2);
    store_rmses_conv(i, :) = rmse_conv;
    
    
    %% Project spline to convex-monotone.
    [y_cm, cm_iter] = convex_monotone_2d(tt, y_nsy, shape);
    y_cm = y_cm(:);

    % Compute RMSE for spline convex, and store results.
    rmse_cm = sqrt(1/length(tt) * norm(y_cm - ytruth_on_test)^2);
    store_rmses_cm(i, :) = rmse_cm;
    
    
    %% Plot results
    if do_plot
        % Plot kernel data.
        figure
        subplot(1, 5, 1)
        mesh(t1, t2, reshape(ytruth_on_test, dim, dim)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zl = zlim;
        title(sprintf('Monotone Truth'))

        subplot(1, 5, 2)
        mesh(t1, t2, reshape(y_nsy, dim, dim)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('Spline Unconstrained (RMSE = %s)', num2str(rmse, '%0.3f')));

        subplot(1, 5, 3)
        mesh(t1, t2, reshape(y_mono, dim, dim)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('Spline Monotone (RMSE = %s)', num2str(rmse_mono, '%0.3f')));
        
        subplot(1, 5, 4)
        mesh(t1, t2, reshape(y_conv, dim, dim)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('Spline Convex (RMSE = %s)', num2str(rmse_conv, '%0.3f')));
        
        subplot(1, 5, 5)
        mesh(t1, t2, reshape(y_cm, dim, dim)); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        zlim(zl);
        title(sprintf('Spline CM (RMSE = %s)', num2str(rmse_cm, '%0.3f')));
    end
    
end

toc(start_time);

%% 
% Save results to file.
results = [store_rmses store_rmses_mono store_rmses_conv store_rmses_cm];
csvwrite(sprintf('data/rmses_spline_orig_mono_conv_cm_50data_%s.csv', shape), results);

if strcmp(platform, 'mac')
    sendmail('momod@utexas.edu', ...
        sprintf('RESULTS mac: shively_%s_rmses', shape), ...
        sprintf('Total time: %s', num2str(toc(start_time), '%0.2f')), ...
        {strcat('/Users/mauricediesendruck/Google Drive/', ...
                '0-LIZHEN RESEARCH/gp/data/', ...
                sprintf('rmses_spline_orig_mono_conv_cm_50data_%s.csv', shape))});       
elseif strcmp(platform, 'linux')
    sendmail('momod@utexas.edu', ...
        sprintf('RESULTS linux: shively_%s_rmses', shape), ...
        sprintf('Total time: %s', num2str(toc(start_time), '%0.2f')), ...
        {strcat('/home/momod/Documents/gp/data/', ...
                sprintf('rmses_spline_orig_mono_conv_cm_50data_%s.csv', shape))});
end

    
    