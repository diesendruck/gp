% Constants.
n = 100;
d = 2;
shape = 'chair';
do_grid = 1;
data_grid_gran = 5;

% Create and show original surface.
[x_nsy, x_nsy_jit, y_nsy, y_nsy_jit] = make_noisy_convex(n, d, ...
shape, do_grid, data_grid_gran);
p3(x_nsy, y_nsy)

% Create and show monotone projection.
