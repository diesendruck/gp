% Test kernel, kernel projection, and Sen methods on sample data sets.

% Constants.
n_data = 100;  % Number of data points.
d = 2;
do_grid = 1;
data_grid_gran = 10;
mesh_gran = data_grid_gran*2;
shape = 'exponential';
eps1 = 1e-4;              
eps2 = 1e-4; 

% Import data.
filename = 'sample_exponential_data.mat';
xydat = importdata(filename);
x_nsy = xydat(:, [1:2]);
y_samples = xydat(:, [3:size(xydat, 2)]);
num_samples = size(y_samples, 2);

% Set up storage variable for each method, reporting RMSEs.
rmse_kern = zeros(num_samples, 1);
rmse_gp = zeros(num_samples, 1);
rmse_kern_m = zeros(num_samples, 1);
rmse_gp_m = zeros(num_samples, 1);


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


%%
% For each set of y outputs, produce the kernel estimate, the kernel
% projection estimate, and the Sen estimate.
for i = 1:num_samples
    % Get y sample for this iteration.
    y_nsy = y_samples(:, i);
    
    % Kernel regression.
    h0 = rand(size(x_nsy', 1), 1)*10;
    h = Opt_Hyp_Gauss_Ker_Reg(h0, x_nsy', y_nsy');
    y_kern_test = zeros(size(t1));
    for ii = 1:size(t1, 1)
        for jj = 1:size(t1, 2)
            xk = [t1(ii, jj); t2(ii, jj)];
            y_kern_test(ii, jj) = gaussian_kern_reg(xk, x_nsy', ...
                                                    y_nsy', h);
        end
    end
    
    % Kernel regression, convex projection.
    y_kern_test_proj = project_to_convex(length(tt), d, tt, y_kern_test(:), eps1, eps2);

    % Sen, convex projection.
    y_sen_test_proj = project_to_convex(length(tt), d, tt, y_nsy, eps1, eps2);

    % Plot data.
    figure
    subplot(1, 3, 1)
    mesh(t1, t2, y_kern_test); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    zl = zlim;
    
    subplot(1, 3, 2)
    mesh(t1, t2, reshape(y_kern_test_proj, data_grid_gran, data_grid_gran)); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    zlim(zl);
    
    subplot(1, 3, 3)
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_sen_test_proj, 'b.', 'MarkerSize', 20);
    grid on;
    zlim(zl);

    % Store results for this y sample.
    results(:, 1, i) = y_nsy;
    results(:, 2, i) = y_kern_test(:);
    results(:, 3, i) = y_kern_test_proj;
    results(:, 4, i) = y_sen_test_proj;
    
end

save('test_shively_samples.mat', 'results');


    
    