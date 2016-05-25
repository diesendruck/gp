% This version was authored by Maurice Diesendruck.

% Adapted from 
% (1) Demo_example.m
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html
%
% (2) ExampleCAP.m, noisyConvexMin.m
% Hannah, L., and Dunson, D. (2011). Bayesian nonparametric multivariate 
%   convex regression. Available at:
%   https://github.com/laurenahannah/convex-function,
%   https://github.com/laurenahannah/mbcr

start_time = tic
platform = 'mac';
%platform = 'linux';

% Import GPstuff, set paths, and setup email params.
setup_directories_and_email(platform);


%% SET CONSTANTS.
format long;
tol_thres = 0;
eps1 = 10^-5;              % These 2 epsilons are used for convergence of the
eps2 = 10^-5;              %   convex projection algorithm.
iter = 0;                  % Counter for iterations.

n = 100;                   % Data sample size. Only relevant if do_grid=0.
d = 2;                     % Dimension of data points. 
do_grid = 1;               % Indicator for whether to generate random data, or grid data.
data_grid_gran = 10;        % Number of points per dimension. 10 means 10x10 for d=2.

ls_factor = 0.5;           % Lengthscale factor (proportion of x-range).
mesh_gran = 20;            % Number of ticks on mesh for plotting.

if 1
    num_posteriors = 20;      % Number of posterior mcmc samples to generate.
    desired = 2;              % Number of posterior mcmc samples to use.
    mbcr_burn = 1;            % Number of burn-in for MBCR estimate.
    mbcr_tot = 2;             % Number of total samples for MBCR estimate.
    num_global_iters = 1;      % Number of MSEs to produce per shape.
end

if 0
    num_posteriors = 500;      % Number of posterior mcmc samples to generate.
    desired = 20;              % Number of posterior mcmc samples to use.
    mbcr_burn = 500;            % Number of burn-in for MBCR estimate.
    mbcr_tot = 1000;             % Number of total samples for MBCR estimate.
    num_global_iters = 15;      % Number of MSEs to produce per shape.
end


%% SAVE MSE RESULTS TO FILE.
fid = fopen('Results_2d/mses.csv', 'wt');
fprintf(fid, 'data_shape,gp,gp_proj,kern,kern_proj,sen,cap,mbcr\n');


%% CONDUCT EXPERIMENT ON EACH SHAPE.
% List of shapes to run.
shapes = {'chair', 'parabolic_cylinder', 'wolverine', 'trough', ...
    'paraboloid', 'hand', 'exponential', 'hannah2'};

% Run experiment for each shape.
for ii = 1:num_global_iters
    for shape = shapes
        shape_start_time = tic;
        
        [gp, gp_proj, kern, kern_proj, sen, cap, mbcr, ...
         gp_proj_time_elapsed, mbcr_time_elapsed] = ...
            gp_compare_mses_shape(tol_thres, eps1, eps2, iter, n, ...
                ls_factor, mesh_gran, num_posteriors, desired, d, ...
                shape{1}, fid, mbcr_burn, mbcr_tot, do_grid, ...
                data_grid_gran, platform);
        
        shape_time_elapsed = toc(shape_start_time);
        total_time_elapsed = toc(start_time);
        
        sendmail('momod@utexas.edu', ...
            sprintf('UPDATE %s: MSE Comparisons', platform), ...
            sprintf(strcat('Global iter: %s\n', ...
                           'Data shape:  %s\n', ...
                           'gp:                %s\n', ...
                           'gp_proj:        %s\n', ... 
                           'kern:             %s\n', ...
                           'kern_proj:     %s\n', ...
                           'sen:              %s\n', ...
                           'cap:              %s\n', ...
                           'mbcr:            %s\n', ...
                           'gp_proj_time:           %s\n', ...
                           'mbcr_time:            %s\n', ...
                           'shape_time:             %s\n', ...
                           'total_time_elapsed: %s\n'), ...
                num2str(ii), shape{1}, num2str(gp, '%0.2f'), ...
                num2str(gp_proj, '%0.2f'), num2str(kern, '%0.2f'), ...
                num2str(kern_proj, '%0.2f'), num2str(sen, '%0.2f'),...
                num2str(cap, '%0.2f'), num2str(mbcr, '%0.2f'), ...
                num2str(gp_proj_time_elapsed, '%0.2f'), ...
                num2str(mbcr_time_elapsed, '%0.2f'), ...
                num2str(shape_time_elapsed, '%0.2f'), ...
                num2str(total_time_elapsed, '%0.2f')));
    end
end


%% CLOSE AND SEND RESULTS FILE.
fclose(fid);
   
if strcmp(platform, 'mac')
    sendmail('momod@utexas.edu', 'RESULTS mac: MSE Comparisons', ...
        sprintf('Total time: %s', num2str(total_time_elapsed, '%0.2f')),...
        {'/Users/mauricediesendruck/Google Drive/0-LIZHEN RESEARCH/gp/Results_2d/mses.csv'});
elseif strcmp(platform, 'linux')
    sendmail('momod@utexas.edu', 'RESULTS linux: MSE Comparisons', ...
        sprintf('Total time: %s', num2str(total_time_elapsed, '%0.2f')),...
        {'/home/momod/Documents/gp/Results_2d/mses.csv'});
end

end_time = toc(start_time)



    
