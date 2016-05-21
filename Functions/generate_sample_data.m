% Produces (x1, x2) grid, and associated y-output samples.

% Constants.
n = 100;
d = 2;
do_grid = 1;
data_grid_gran = 10;
shape = 'hand';

% Set up storage.
results = zeros(n, 12);

% Produce x_nsy, which is a fixed grid.
[x_nsy, ~, y_nsy, ~] = make_noisy_convex(n, d, shape, do_grid, data_grid_gran);
results(:, 1) = x_nsy(:, 1);
results(:, 2) = x_nsy(:, 2);

% First of 10 samples of output.
results(:, 3) = y_nsy;

% Remaining 9 samples of output.
for i = 1:9
    [~, ~, y_nsy_i, ~] = make_noisy_convex(n, d, shape, do_grid, data_grid_gran);
    results(:, 3+i) = y_nsy_i;
end

save('sample_hand_data.mat', 'results');
save('sample_hand_data.txt', 'results', '-ascii');
