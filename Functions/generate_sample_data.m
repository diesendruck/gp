% Produces (x1, x2) grid, and associated y-output samples.

% Constants.
n = 100;
d = 2;
do_grid = 1;
data_grid_gran = 10;
n_samples = 50;
shape = 'exponential';  % CHOOSE DESIRED SHAPE.

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

save(sprintf('data/shively_%s_samples.dat', shape), 'results');

