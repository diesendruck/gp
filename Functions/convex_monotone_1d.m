function [y_cm, r_mono, r_conv] = convex_monotone_1d(x_nsy, y_nsy)
% Projects 1d function to convex monotone (cm).
%
% Args:
%   x_nsy: n x 1 matrix of data values.
%   y_nsy: n x 1 matrix of response values.
%
% Returns:
%   y_cm: n x 1 matrix of convex monotone values for x_nsy inputs.

% Define convergence indicator and max number of iterations in while loop.
iter = 0;
max_iter = 1e3;
eps = 1e-6;

% Set up initial monotone and convex projections, and their residuals.
f_init = y_nsy;

f_mono = monotone(f_init);
r_mono = f_mono - f_init;

f_conv = convex1(x_nsy, f_init + r_mono);
r_conv = f_conv - (f_init + r_mono);

while (iter < max_iter)
    % If latest convex and monotone functions are the same (epsilon-close),
    % we have the desired convex monotone projection.
    latest_f_m = f_mono(:, max(size(f_mono, 2)));
    latest_f_c = f_conv(:, max(size(f_conv, 2)));
    if 1/length(x_nsy) * norm(latest_f_c - latest_f_m)^2 < eps
       break
    end
    
    % Get previous convex residual.
    latest_r_c = r_conv(:, max(size(r_conv, 2)));
    
    % Compute monotone projection, and its residual.
    new_f_m = monotone_1d(f_init + latest_r_c);
    new_r_m = new_f_m - (f_init + latest_r_c);
    
    % Compute convex projection (of above monotone), and its residual.
    new_f_c = convex_1d(x_nsy, f_init + new_r_m);
    new_r_c = new_f_c - (f_init + new_r_m);
    
    % Store each monotone and convex function, and associated residuals.
    f_mono = [f_mono new_f_m];
    r_mono = [r_mono new_r_m];
    f_conv = [f_conv new_f_c];
    r_conv = [r_conv new_r_c];
    
    iter = iter + 1;
end

% Return one of the functions.
y_cm = latest_f_m;

% Plot original data and cm projection.
plot(x_nsy, y_nsy, 'kx', 'markersize', 10); hold on
plot(x_nsy, latest_f_m, 'r-', 'markersize', 10)
plot(x_nsy, latest_f_m, 'r.', 'markersize', 10)

sprintf('Used %d iterations', iter)

end
