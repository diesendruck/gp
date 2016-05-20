% Convex and Monotone Projections of Gaussian Processes.

start_time = tic;

%% CONSTANTS
tol_thres = 0;
eps1 = 1e-5;                % These 2 epsilons are used for convergence of the
eps2 = 1e-5;                %   convex projection algorithm.
iter = 0;                   % Counter for iterations.

n = 100;                    % Data sample size. Only relevant if do_grid=0.
d = 2;                      % Dimension of data points. 
do_grid = 1;                % Indicator for whether to generate random data, or grid data.
data_grid_gran = 10;        % Number of points per dimension. 10 means 10x10 for d=2.
dim = data_grid_gran;
mesh_gran = 2*dim;          % Number of ticks on mesh for plotting.
ls_factor = 0.5;            % Lengthscale factor (proportion of x-range).

verbose = 1;                % Whether to print things to console.
do_plot = 1;                % Whether to make plots, versus just results by email.

% Short run.
if 1
    num_posteriors = 50;    % Number of posterior mcmc samples to generate.
    desired = 2;            % Number of posterior mcmc samples to use.
    num_global_iters = 1;   % Number of MSEs to produce per shape.
end

% Full run.
if 0
    num_posteriors = 2000;  % Number of posterior mcmc samples to generate.
    desired = 50;           % Number of posterior mcmc samples to use.
    num_global_iters = 15;  % Number of MSEs to produce per shape.
end

% Choose platform 'mac' or 'linux'.
platform = 'mac';

% Choose 2D shape.
shape = 'cm1';


%% IMPORT GPSTUFF AND SET PATHS.
if strcmp(platform, 'mac');  % Mac version.
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/GPstuff-4.6/
    matlab_install
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/Programs/
    run_mex_commands();
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Programs')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Functions')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Smoothing')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/convex-function')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/mbcr')
elseif strcmp(platform, 'linux');  % Linux version.
    cd ~/Documents/gp/GPstuff-4.6/
    matlab_install
    cd ~/Documents/gp/Programs/
    run_mex_commands();
    cd ~/Documents/gp/
    addpath('~/Documents/gp/')
    addpath('~/Documents/gp/Programs')
    addpath('~/Documents/gp/Functions')
    addpath('~/Documents/gp/Smoothing')
    addpath('~/Documents/gp/convex-function')
    addpath('~/Documents/gp/mbcr')
end


%% EMAIL PARAMS
myaddress = 'eltegog@gmail.com';
myp = 'T0g.eltegog';
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',myp);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');


%% SET UP FILE FOR RESULTS
fid = fopen('Results_2d/cm_mses.csv', 'wt');
fprintf(fid, 'data_shape,gp,gp_conv,gp_mono,gp_cm,kern,kern_conv,kern_mono,kern_cm\n');

    
%% GLOBAL RUN
for i = 1:num_global_iters;   
    %% MAKE TEST DATA
    [x_nsy, ~, y_nsy, ~] = make_noisy_convex(n, d, shape, do_grid, ...
        data_grid_gran);


    %% ORIGINAL CONVEX MONOTONE DATA
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

    if do_plot
        figure
        x0=50;
        y0=100;
        width=1200;
        height=600;
        set(gcf,'units','points','position',[x0,y0,width,height])

        %figure('units','normlized','position',[.1 .1 .4 .4])

        yq_conv = griddata(xt(:, 1), xt(:, 2), ytruth_on_grid, xt1, xt2);
        subplot(2, 5, 1);
        surf(xt1, xt2, yq_conv); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        title('True Convex');

        % Store original z-axis limits.
        zl = zlim;
    end


    %% GP POSTERIORS, THEIR CONVEX PROJECTIONS, AVERAGES, AND MSESs
    gp_c_time = tic;

    % Get samples from GP posterior MCMC, project each to convex, and store.
    [xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy, ...
        y_nsy, ls_factor, num_posteriors, mesh_gran);
    n_gp = length(xt);

    n_entries = min(desired, posterior_sample_count); 
    mcmcs = zeros(n_gp, n_entries);
    convs = zeros(n_gp, n_entries);

    for index = 1:n_entries
        if verbose
            fprintf('(%s) Projecting mcmc surface to convex, sample %d.\n', ...
                shape, index);
        end
        % Sample once from posterior, and store it as a column in mcmcs.
        y_smp = Eft_s(:, randi(posterior_sample_count));
        mcmcs(:, index) = y_smp;
        % Get convex projection of sample, and store it as a column in projs.
        c_time = tic;
        y_smp_convex = project_to_convex(n_gp, d, xt, y_smp, eps1, eps2);
        fprintf('(Convex 2D) Time: %d\n', toc(c_time))
        convs(:, index) = y_smp_convex;
    end

    % Compute averages over mcmc and convex projections, respectively.
    avg_mcmcs = mean(mcmcs, 2);  % Row means.
    avg_convs = mean(convs, 2);

    % Evaluate mcmc and proj estimates on test points.
    y_gp_test = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, t1, t2);
    y_gp_conv_test = griddata(xt(:, 1), xt(:, 2), avg_convs, t1, t2);
    mse_gp = 1/length(tt) * norm(y_gp_test(:) - ytruth_on_test)^2;
    mse_gp_conv = 1/length(tt) * norm(y_gp_conv_test(:) - ytruth_on_test)^2;

    if do_plot
        % Plot avg MCMC over original data.
        subplot(2, 5, 2);
        yq_mcmc = griddata(xt(:, 1), xt(:, 2), avg_mcmcs, xt1, xt2);
        surf(xt1, xt2, yq_mcmc); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        title(sprintf('AvgGP (MSE = %s)', num2str(mse_gp, '%0.3f')));

        % Set z-axis limits to match original plot.
        zlim(zl);

        % Plot avg convex proj over original data.
        subplot(2, 5, 3); 
        yq_conv = griddata(xt(:, 1), xt(:, 2), avg_convs, xt1, xt2);
        surf(xt1, xt2, yq_conv); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        title(sprintf('AvgGP Convex (MSE = %s)', num2str(mse_gp_conv, '%0.3f')));

        % Set z-axis limits to match original plot.
        zlim(zl);
    end

    % Report time of subroutine.
    fprintf('(AvgGP and AvgGP Convex) Time: %d\n', toc(gp_c_time))


    %% MONOTONE PROJECTION OF AVG GP.
    % Time subroutine.
    gp_m_time = tic;

    % Do monotone projection of Avg GP, on full mesh "xt".
    %f = monotone_2d(x_nsy, y_mcmc_test);
    %f_mono = griddata(tt(:, 1), tt(:, 2), f(:), xt1, xt2);
    y_gp_mono = monotone_2d(xt, reshape(avg_mcmcs, mesh_gran, mesh_gran));

    % Compute MSE only over test points.
    y_gp_mono_test = griddata(xt(:, 1), xt(:, 2), y_gp_mono(:), t1, t2);
    mse_gp_mono = 1/length(tt) * norm(y_gp_mono_test(:) - ytruth_on_test)^2;

    subplot(2, 5, 4);
    surf(xt1, xt2, y_gp_mono); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('AvgGP Monotone (MSE = %s)', num2str(mse_gp_mono, '%0.3f')));

    % Set z-axis limits to match original plot.
    zlim(zl);

    % Report time of subroutine.
    fprintf('(AvgGP Monotone 2D) Time: %d\n', toc(gp_m_time))


    %% CONVEX MONOTONE PROJECTION OF AVG GP.
    % Time subroutine.
    gp_cm_time = tic;

    % Do convex monotone projection of Avg GP, on full mesh "xt".
    %f = convex_monotone_2d(x_nsy, y_mcmc_test);
    %f_cm = griddata(tt(:, 1), tt(:, 2), f(:), xt1, xt2);
    y_gp_cm = convex_monotone_2d(xt, reshape(avg_mcmcs, mesh_gran, mesh_gran));

    % Compute MSE only over test points.
    y_gp_cm_test = griddata(xt(:, 1), xt(:, 2), y_gp_cm(:), t1, t2);
    mse_gp_cm = 1/length(tt) * norm(y_gp_cm_test(:) - ytruth_on_test)^2;

    subplot(2, 5, 5);
    surf(xt1, xt2, y_gp_cm); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('AvgGP Conv+Mono (MSE = %s)', num2str(mse_gp_cm, '%0.3f')));

    % Set axis limits to match original plot.
    xlim(xl); ylim(yl); zlim(zl);

    % Report time of subroutine.
    fprintf('(AvgGP Convex Monotone 2D) Time: %d\n', toc(gp_cm_time))


    %% KERNEL REGRESSION, ITS CONVEX PROJECTION, AND MSES.
    % Time subroutine.
    kern_c_time = tic;

    % --------FULL MESH--------
    % Kernel regression on full grid mesh "xt".
    h0 = rand(size(x_nsy', 1), 1)*10;
    if verbose
        fprintf('(%s) Optimizing for kernel regression bandwidth.\n', shape);
    end
    h = Opt_Hyp_Gauss_Ker_Reg(h0, x_nsy', y_nsy');
    y_kern = zeros(size(xt1));
    for ii = 1:size(xt1, 1)
        for jj = 1:size(xt1, 2)
            xk = [xt1(ii, jj); xt2(ii, jj)];
            y_kern(ii, jj) = gaussian_kern_reg(xk, x_nsy', y_nsy', h);
        end
    end

    % Do convex projection of kernel regression.
    y_kern_conv = project_to_convex(length(xt), d, xt, y_kern(:), eps1, eps2);

    % Evaluate kern and kern_conv estimates on test points.
    y_kern_test = griddata(xt(:, 1), xt(:, 2), y_kern(:), t1, t2);
    y_kern_conv_test = griddata(xt(:, 1), xt(:, 2), y_kern_conv(:), t1, t2);

    % Compute mses on test data.
    mse_kern = 1/length(tt) * norm(y_kern_test(:) - ytruth_on_test)^2;
    mse_kern_conv = 1/length(tt) * norm(y_kern_conv_test(:) - ytruth_on_test)^2;


    % Plot kernel regression and convex projection over original data.
    if do_plot
        subplot(2, 5, 7);
        yq_kern = griddata(xt(:, 1), xt(:, 2), y_kern(:), xt1, xt2);
        surf(xt1, xt2, yq_kern); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        title(sprintf('Kernel (MSE = %s)', num2str(mse_kern, '%0.3f')));
        % Standardize z-axis.
        zlim(zl);

        subplot(2, 5, 8);
        yq_kern_conv = griddata(xt(:, 1), xt(:, 2), y_kern_conv(:), xt1, xt2);
        surf(xt1, xt2, yq_kern_conv); hold on;
        plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
        title(sprintf('Kernel Convex (MSE = %s)', num2str(mse_kern_conv, '%0.3f')));
        % Standardize z-axis.
        zlim(zl);
    end

    % Report time of subroutine.
    fprintf('(Kernel and Kernel Convex 2D) Time: %d\n', toc(kern_c_time))


    %% MONOTONE PROJECTION OF KERNEL REGRESSION.
    % Time subroutine.
    kern_m_time = tic;

    % Do monotone projection of kernel regression, on full mesh "xt".
    y_kern_mono = monotone_2d(xt, y_kern);

    % Compute MSE only over test points.
    y_kern_mono_test = griddata(xt(:, 1), xt(:, 2), y_kern_mono(:), t1, t2);
    mse_kern_mono = 1/length(tt) * norm(y_kern_mono_test(:) - ytruth_on_test)^2;

    subplot(2, 5, 9);
    surf(xt1, xt2, y_kern_mono); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('Kernel Monotone (MSE = %s)', num2str(mse_kern_mono, '%0.3f')));

    % Set z-axis limits to match original plot.
    zlim(zl);

    % Report time of subroutine.
    fprintf('(Kernel Monotone 2D) Time: %d\n', toc(kern_m_time))


    %% CONVEX MONOTONE PROJECTION OF KERNEL REGRESSION.
    % Time subroutine.
    kern_cm_time = tic;

    % Do convex monotone projection of Avg GP, on full mesh "xt".
    %f = convex_monotone_2d(x_nsy, y_mcmc_test);
    %f_cm = griddata(tt(:, 1), tt(:, 2), f(:), xt1, xt2);
    y_kern_cm = convex_monotone_2d(xt, y_kern);

    % Compute MSE only over test points.
    y_kern_cm_test = griddata(xt(:, 1), xt(:, 2), y_kern_cm(:), t1, t2);
    mse_kern_cm = 1/length(tt) * norm(y_kern_cm_test(:) - ytruth_on_test)^2;

    subplot(2, 5, 10);
    surf(xt1, xt2, y_kern_cm); hold on;
    plot3(x_nsy(:, 1), x_nsy(:, 2), y_nsy, 'r.', 'MarkerSize', 20);
    title(sprintf('Kernel Conv+Mono (MSE = %s)', num2str(mse_kern_cm, '%0.3f')));

    % Set axis limits to match original plot.
    xlim(xl); ylim(yl); zlim(zl);

    % Report time of subroutine.
    fprintf('(Kernel Convex Monotone 2D) Time: %d\n', toc(kern_cm_time))

    
    %% ADD SUMMARY TEXT TO PLOT.
    % Add text on plot to say which method did best (lowest MSE).
    if do_plot
        ax = subplot(2, 5, 6);
        [~, index] = min([mse_gp mse_gp_conv mse_gp_mono mse_gp_cm ...
                          mse_kern mse_kern_conv mse_kern_mono mse_kern_cm]);
        methods = {'mse_gp' 'mse_gp_conv' 'mse_gp_mono' 'mse_gp_cm' ...
                   'mse_kern' 'mse_kern_conv' 'mse_kern_mono' 'mse_kern_cm'};
        min_str = strrep(char(methods(index)), '_', '\_');
        text(0, 0.5, 'Min MSE Method:', 'FontSize', 14);
        text(0, 0.3, min_str, 'FontSize', 14);
        set (ax, 'visible', 'off')

    end

    
    %% SAVE FILE DATA AND FIGURE.
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s\n', shape, ...
        num2str(mse_gp, '%0.7f'), ...
        num2str(mse_gp_conv, '%0.7f'), ...
        num2str(mse_gp_mono, '%0.7f'), ...
        num2str(mse_gp_cm, '%0.7f'), ...
        num2str(mse_kern, '%0.7f'), ...
        num2str(mse_kern_conv, '%0.7f'), ...
        num2str(mse_kern_mono, '%0.7f'), ...
        num2str(mse_kern_cm, '%0.7f')); 
    
end


%% CLOSE AND SEND RESULTS FILE.
fclose(fid);
   
if strcmp(platform, 'mac')
    sendmail('momod@utexas.edu', 'RESULTS mac: CM_MSE', ...
        sprintf('Total time: %s', num2str(toc(start_time), '%0.2f')), ...
        {'/Users/mauricediesendruck/Google Drive/0-LIZHEN RESEARCH/gp/Results_2d/cm_mses.csv'});
elseif strcmp(platform, 'linux')
    sendmail('momod@utexas.edu', 'RESULTS linux: CM_MSE', ...
        sprintf('Total time: %s', num2str(toc(start_time), '%0.2f')), ...
        {'/home/momod/Documents/gp/Results_2d/cm_mses.csv'});
end


%%
toc(start_time)
