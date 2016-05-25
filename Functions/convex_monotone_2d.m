function [y_cm, iter] = convex_monotone_2d(x_nsy, y_nsy, shape)
% Projects 2D data to convex monotone (cm).
%
% Args:
%   x_nsy: Matrix (nx2) of data input values.
%   y_nsy: Matrix (nx1) of noisy response values.
%   shape: String, e.g. "cm1" or "hand".
%
% Returns:
%   y_cm: Matrix (sqrt(n)xsqrt(n)) of convex monotone response variable.
%   iter: Integer of the iteration at which the algorithm converged.


%% SETUP

% Toggle plotting.
do_plot = 0;

% Define convergence indicator and max number of iterations of while loop.
iter = 1;
max_iter = 1e3;
eps = 1e-6;

% Set up initial surface.
dim = sqrt(length(x_nsy));
x1 = reshape(x_nsy(:, 1), dim, dim)';
x2 = reshape(x_nsy(:, 2), dim, dim)';
xt = [x1(:) x2(:)];
ytruth = compute_truth_from_xt(xt, shape);
f_init = reshape(y_nsy, dim, dim)';

% Show original surface.
if do_plot
    x0=50;
    y0=100;
    width=1200;
    height=400;
    set(gcf,'units','points','position',[x0,y0,width,height])

    subplot(2, 4, 1);
    surf(x1, x2, f_init); hold on; 
    p3(x_nsy, ytruth);
    title('Initial Surface');
end


%% FIRST ITERATION OF PROJECTIONS

% First monotone row projection (of initial matrix), and its residual.
f_row = zeros(size(f_init));
for row = 1:size(f_init, 1)
    f_row(row, :) = monotoneUp_1d(f_init(row, :));
end
r_row = f_row - f_init;

% First monotone column projection (of monotone row projection above), and
% its residual.
f_col = zeros(size(f_init));
f = f_init + r_row;  % Initial, plus row residual.
for col = 1:size(f_init, 2)
    f_col(:, col) = monotoneUp_1d(f(:, col));
end
r_col = f_col - f;

% First convex projection (of monotone column projection above), and its 
% residual.
f = f_init + r_col + r_row;  % Initial, plus row and col residuals.
f_conv = convex_2d(x_nsy, f(:));
f_conv = reshape(f_conv, dim, dim);
r_conv = f_conv - f;

% Compute MSE after first set of projections.
mse_iter_1 = 1/length(xt) * norm(f_conv(:) - ytruth)^2;

% Show surface after first row and column projection.
if do_plot
    subplot(2, 4, 2);
    surf(x1, x2, f_row); hold on;
    p3(x_nsy, ytruth);
    title(sprintf('Row proj, iter: %s', num2str(iter)));

    subplot(2, 4, 3);
    surf(x1, x2, f_col); hold on;
    p3(x_nsy, ytruth);
    title(sprintf('Column proj, iter: %s', num2str(iter)));

    subplot(2, 4, 4);
    surf(x1, x2, f_conv); hold on;
    p3(x_nsy, ytruth);
    title(sprintf('Convex, iter: %s. MSE=%s', num2str(iter), ...
                  num2str(mse_iter_1, '%0.2f')));
end


%% OPTIMIZATION LOOP

% Alternately do row and col projection until convergence.
while (iter < max_iter)
    % If latest row and column matrices are the same (epsilon-close), we 
    % have the desired 2D monotone projection.
    latest_f_row = f_row(:, :, size(f_row, 3));
    latest_f_col = f_col(:, :, size(f_col, 3));
    if 1/numel(f_init) * norm(latest_f_row - latest_f_col)^2 < eps
        disp(sprintf('(Convex Monotone 2D) Converged at iter %s', ...
            num2str(iter)))
        break
    end
    
    % Get previous residuals.
    latest_r_col = r_col(:, :, size(r_col, 3));
    latest_r_conv = r_conv(:, :, size(r_conv, 3));
    
    % Get new ROW PROJECTION, and its residual.
    % Pre-allocate space for proj.
    new_f_row = zeros(size(f_init));
    % Use initial, plus latest column and convex residuals.
    f = f_init + latest_r_col + latest_r_conv;
    for row = 1:size(f_init, 1)
        new_f_row(row, :) = monotoneUp_1d(f(row, :));
    end
    new_r_row = new_f_row - f;
    
    % Get new COLUMN PROJECTION, and its residual.
    % Pre-allocate space for proj.
    new_f_col = zeros(size(f_init));
    % Use initial, plus latest convex residual and row residual above.
    f = f_init + latest_r_conv + new_r_row;
    for col = 1:size(f_init, 2)
        new_f_col(:, col) = monotoneUp_1d(f(:, col));
    end
    new_r_col = new_f_col - f;
    
    % Get new CONVEX PROJECTION, and its residual.
    % Use initial, plus row residual above and column residual above.
    f = f_init + new_r_row + new_r_col;
    new_f_conv = convex_2d(x_nsy, f(:));
    new_f_conv = reshape(new_f_conv, dim, dim);
    new_r_conv = new_f_conv - f;
    
    % Store each row and col projection, and associated residuals.
    f_row(:, :, iter+1) = new_f_row;
    r_row(:, :, iter+1) = new_r_row;
    f_col(:, :, iter+1) = new_f_col;
    r_col(:, :, iter+1) = new_r_col;
    f_conv(:, :, iter+1) = new_f_conv;
    r_conv(:, :, iter+1) = new_r_conv;
    
    iter = iter + 1;
    
end

% Compute MSE for convex projection, after projections converge.
y_cm = f_conv(:, :, size(f_conv, 3));
mse_iter_converge = 1/length(xt) * norm(y_cm(:) - ytruth)^2;


%% SHOW RESULT

% Show last projection.
if do_plot
    subplot(2, 4, 6);
    f = f_row(:, :, size(f_row, 3));
    surf(x1, x2, f); hold on;
    p3(x_nsy, ytruth);
    title(sprintf('Row proj, iter: %s', num2str(iter)));

    subplot(2, 4, 7);
    f = f_col(:, :, size(f_col, 3));
    surf(x1, x2, f); hold on;
    p3(x_nsy, ytruth);
    title(sprintf('Column proj, iter: %s', num2str(iter)));

    subplot(2, 4, 8);
    f = f_conv(:, :, size(f_conv, 3));
    surf(x1, x2, f); hold on;
    p3(x_nsy, ytruth);
    title(sprintf('Convex, iter: %s. MSE=%s', num2str(iter), ...
                      num2str(mse_iter_converge, '%0.2f')));
end

end
