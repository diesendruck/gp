% Script to produce all plots in one.

n = 100;
d = 1;
do_grid = 1;
data_grid_gran = 100;
mesh_gran = 200;

shapes = {'1d_parabola', '1d_exponential', '1d_negative_entropy'};

figure; hold on;

for i = 1:length(shapes)
    shape = shapes{i};
    
    % Setup storage for y values and function values. 
    if do_grid
        len = data_grid_gran;
    else
        len = n;
    end
    f = zeros(len, 1);

    % Define x values.
    if strcmp(shape, '1d_parabola');
        if do_grid
            x_nsy = linspace(0, 10, len)';
        else
            x_nsy = unifrnd(0, 10, len, 1);
        end

        sig = 2.0;
        noise = sig*randn(len, 1);

        % Compute function value.
        for ii = 1:len,
            f(ii) = 0.004 * (x_nsy(ii)-7)^4;
        end

    elseif strcmp(shape, '1d_exponential');
        if do_grid
            x_nsy = linspace(0, 10, len)';
        else
            x_nsy = unifrnd(0, 10, len, 1);
        end

        sig = 2.0;
        noise = sig*randn(len, 1);

        % Compute function value.
        for ii = 1:len,
            %f(ii) = exp(x_nsy(ii)) + 0.5*x_nsy(ii)^2;
            f(ii) = 4.5e-4 * exp(x_nsy(ii));
        end

    elseif strcmp(shape, '1d_negative_entropy');
        if do_grid
            x_nsy = linspace(0, 10, len)';
        else
            x_nsy = unifrnd(0, 10, len, 1);
        end

        sig = 2.0;
        noise = sig*randn(len,1);

        % Compute function value.
        for ii = 1:len,
            f(ii) = 0.4 * x_nsy(ii) * log(x_nsy(ii)+1);
        end

    else
        error('Shape not recognized.')
    end

    % Response values = convex + noise.
    y_nsy = f + noise;

    % Plot the true surface, with noisy points on top.
    [~, ~, ~, x_grid] = compute_mesh_info_1d(x_nsy, ...
        mesh_gran);
    subplot(1, 3, i);
    yq = interp1(x_nsy, f, x_grid);
    plot(x_grid, yq); hold on;
    plot(x_nsy, y_nsy, 'r.', 'MarkerSize', 10);
    title(strrep(shape, '_', '\_')); hold on;
    ylim([-6 14])
end
