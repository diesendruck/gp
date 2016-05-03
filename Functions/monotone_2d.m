function Pmono = monotone_2d(x_nsy, y_nsy)
% Project vector to monotone space.
%
% Args:
%   x_nsy: Matrix (nx2) of data input values.
%   y_nsy: Matrix (nx1) of noisy response values.
%
% Returns:
%   Pmono: Matrix (sqrt(n)xsqrt(n)) of the monotone response variable.

% Toggle plotting.
do_plot = 1;

% Define convergence indicator and max number of iterations of while loop.
iter = 1;
max_iter = 1e2;
eps = 1e-6;

% Set up initial surface.
dim = sqrt(length(x_nsy));
x1 = reshape(x_nsy(:, 1), dim, dim)';
x2 = reshape(x_nsy(:, 2), dim, dim)';
xt = [x1(:) x2(:)];
ytruth = compute_truth_from_xt(xt, 'cm1');
f_init = reshape(y_nsy, dim, dim)';

% Show original surface.
if do_plot
    figure; surf(x1, x2, f_init); hold on; 
    p3(x_nsy, ytruth);
    %scatter3(x_nsy(:,1), x_nsy(:,2), ytruth, ones(100, 1)*30, 'r', 'filled')
end

% First monotone row projection, and its residual.
f_row = zeros(size(f_init));
for row = 1:size(f_init, 1)
    f_row(row, :) = monotoneUp_1d(f_init(row, :));
end
r_row = f_init - f_row;

% First monotone column projection, and its residual.
f_col = zeros(size(f_init));
f = f_init + r_row;  % Add back row residual.
for col = 1:size(f_init, 2)
    f_col(:, col) = monotoneUp_1d(f(:, col));
end
r_col = f_row - f_col;

% Show surface after first row and column projection.
if do_plot
    figure; surf(x1, x2, f_col); hold on;
    p3(x_nsy, ytruth);
end

% Alternately do row and col projection until convergence.
while (iter < max_iter)
    % If latest row and column matrices are the same (epsilon-close), we 
    % have the desired 2D monotone projection.
    latest_f_row = f_row(:, :, max(size(f_row, 3)));
    latest_f_col = f_col(:, :, max(size(f_col, 3)));
    if 1/numel(f_init) * norm(latest_f_row - latest_f_col)^2 < eps
        sprintf('Converged at iter %s', num2str(iter))
        break
    end
    
    % Get previous column residual.
    latest_r_col = r_col(:, :, max(size(r_col, 3)));
    
    % Get row projection, and its residual.
    new_f_row = zeros(size(f_init));
    f = f_init + latest_r_col;
    for row = 1:size(f_init, 1)
        new_f_row(row, :) = monotoneUp_1d(f(row, :));
    end
    new_r_row = new_f_row - f;
    
    % Get col projection (with row residual from above), and its residual.
    new_f_col = zeros(size(f_init));
    f = f_init + new_r_row;
    for col = 1:size(f_init, 2)
        new_f_col(:, col) = monotoneUp_1d(f(:, col));
    end
    new_r_col = new_f_col - f;
    
    % Store each row and col projection, and associated residuals.
    f_row(:, :, iter+1) = new_f_row;
    r_row(:, :, iter+1) = new_r_row;
    f_col(:, :, iter+1) = new_f_col;
    r_col(:, :, iter+1) = new_r_col;
    
    if do_plot,
        if mod(iter,10)==0
            figure; surf(x1, x2, f_col(:, :, iter)); hold on;
            p3(x_nsy, ytruth);
        end
    end
    
    iter = iter + 1;
    
end

% Show last projection.
if do_plot
    f = f_col(:, :, max(size(f_col, 3)));
    figure; surf(x1, x2, f); hold on;
    p3(x_nsy, ytruth);
    %scatter3(x_nsy(:,1), x_nsy(:,2), ytruth, ones(100, 1)*30, 'r', 'filled')
end

end