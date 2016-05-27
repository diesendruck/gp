% Produces (x1, x2) grid, and associated y-output samples.
% Files stored in gp/data/ with the format: 
% 50data_samples_[shape name].csv

% Constants.
n = 100;
d = 2;
do_grid = 1;
data_grid_gran = 10;
n_samples = 50;
shapes = {'chair', 'parabolic_cylinder', 'wolverine', 'trough', ...
    'paraboloid', 'hand', 'exponential', 'hannah2'};
shapes = {'cm1', 'cm2', 'cm3', 'cm4'};

for shape = shapes
    shape = shape{1};
    % Set up storage.
    results = zeros(n, n_samples+2);

    % Produce x_nsy, which is a fixed grid.
    [x_nsy, ~, y_nsy, ~] = make_noisy_convex(n, d, shape, do_grid, data_grid_gran);
    results(:, 1) = x_nsy(:, 1);
    results(:, 2) = x_nsy(:, 2);

    % First of n_samples outputs.
    results(:, 3) = y_nsy;

    % Remaining samples of output.
    for i = 1:n_samples-1
        [~, ~, y_nsy_i, ~] = make_noisy_convex(n, d, shape, do_grid, data_grid_gran);
        results(:, 3+i) = y_nsy_i;
    end

    csvwrite(sprintf('data/50data_samples_%s.csv', shape), results);
    
end
