% This version was partially authored by Maurice Diesendruck.

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


%% SETUP EMAIL PARAMS.
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


%% SET CONSTANTS.
format long;
tol_thres = 0;
eps1 = 10^-5;              % These 2 epsilons are used for convergence of the
eps2 = 10^-5;              %   convex projection algorithm.
iter = 0;                  % Counter for iterations.

n = 100;                   % Data sample size. Ignored if do_grid=1.
d = 1;                     % Dimension of data points.
do_grid = 1;               % Indicator for whether to generate random data, or grid data.
data_grid_gran = 20;       % Number of points per dimension. 10 means 10x1 for d=1.

gp_optimization = 'mcmc';  % "mcmc" or "map" to select GP optimization type. 
                           % TODO: MAP not adjusted for test data.
ls_factor = 0.8;           % Lengthscale factor ([0, 1] proportion of x-range).
mesh_gran = 40;            % Number of ticks on mesh for plotting.

if 1
    num_posteriors = 120;      % Number of posterior mcmc samples to generate.
    desired = 5;              % Number of posterior mcmc samples to use.
    mbcr_burn = 5;            % Number of burn-in for MBCR estimate.
    mbcr_tot = 10;             % Number of total samples for MBCR estimate.
    num_global_iters = 1;      % Number of MSEs to produce per shape.
end

if 0
    num_posteriors = 2000;      % Number of posterior mcmc samples to generate.
    desired = 50;              % Number of posterior mcmc samples to use.
    mbcr_burn = 500;            % Number of burn-in for MBCR estimate.
    mbcr_tot = 1000;             % Number of total samples for MBCR estimate.
    num_global_iters = 5;      % Number of MSEs to produce per shape.
end


%% SAVE MSE RESULTS TO FILE.
fid = fopen('Results_1d/mses.csv', 'wt');
fprintf(fid, 'data_shape,gp,gp_proj,kern,kern_proj,sen,cap,mbcr\n');

%% CONDUCT EXPERIMENT ON EACH SHAPE.
% List of shapes to run.
shapes = {'parabola', 'exponential', 'negative_entropy'};

% Run whole experiment many times.
for ii = 1:num_global_iters
    for shape = shapes
        shape_start_time = tic;

        [gp, gp_proj, kern, kern_proj, sen, cap, mbcr] = ...
            gp_compare_mses_shape_1d(tol_thres, eps1, eps2, iter, n, ...
                ls_factor, mesh_gran, gp_optimization, num_posteriors, ...
                desired, d, shape{1}, fid, mbcr_burn, mbcr_tot, do_grid,...
                data_grid_gran);
        
        shape_time_elapsed = toc(shape_start_time);
        total_time_elapsed = toc(start_time);
                
        sendmail('momod@utexas.edu', ...
            sprintf('UPDATES_1D %s: MSE Comparisons', platform), ...
            sprintf(strcat('Global iter: %s\n', ...
                           'Data shape:  %s\n', ...
                           'gp:                %s\n', ...
                           'gp_proj:        %s\n', ... 
                           'kern:             %s\n', ...
                           'kern_proj:     %s\n', ...
                           'sen:              %s\n', ...
                           'cap:              %s\n', ...
                           'mbcr:            %s\n', ...
                           'shape_time:             %s\n', ...
                           'total_time_elapsed: %s\n'), ...
                num2str(ii), shape{1}, num2str(gp, '%0.2f'), ...
                num2str(gp_proj, '%0.2f'), num2str(kern, '%0.2f'), ...
                num2str(kern_proj, '%0.2f'), num2str(sen, '%0.2f'),...
                num2str(cap, '%0.2f'), ...
                num2str(mbcr, '%0.2f'), ...
                num2str(shape_time_elapsed, '%0.2f'), ...
                num2str(total_time_elapsed, '%0.2f')));
    end
end


%% CLOSE AND SEND RESULTS FILE.
fclose(fid);
   
if strcmp(platform, 'mac')
    sendmail('momod@utexas.edu', 'RESULTS_1D mac', ...
        sprintf('Total time: %s', num2str(total_time_elapsed, '%0.2f')),...
        {'/Users/mauricediesendruck/Google Drive/0-LIZHEN RESEARCH/gp/Results_1d/mses.csv'});
elseif strcmp(platform, 'linux')
    sendmail('momod@utexas.edu', 'RESULTS_1D linux', ...
        sprintf('Total time: %s', num2str(total_time_elapsed, '%0.2f')),...
        {'/home/momod/Documents/gp/Results_1d/mses.csv'});
end

end_time = toc(start_time)




