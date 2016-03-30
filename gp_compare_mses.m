% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html

start_time = tic

%% IMPORT GPSTUFF AND SET PATHS.
if 0  % Mac version.
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/GPstuff-4.6/
    matlab_install
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/Programs/
    run_mex_commands();
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Programs')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Functions')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Smoothing')
end

if 0  % Linux version.
    cd ~/Documents/gp/GPstuff-4.6/
    matlab_install
    cd ~/Documents/gp/Programs/
    run_mex_commands();
    cd ~/Documents/gp/
    addpath('~/Documents/gp/')
    addpath('~/Documents/gp/Programs')
    addpath('~/Documents/gp/Functions')
    addpath('~/Documents/gp/Smoothing')
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
eps1 = 10^-5;          % These 2 epsilons are used for convergence of the
eps2 = 10^-5;          %   convex projection algorithm.
iter = 0;              % Counter for iterations.
n = 40;                % Data sample size.
d = 2;                  % Dimension of data points.
ls_factor = 0.01;      % Lengthscale factor (proportion of x-range).
mesh_gran = 10;        % Number of ticks on mesh for plotting.
num_posteriors = 120; % Number of posterior mcmc samples to generate.
desired = 2;         % Number of posterior mcmc samples to use.
num_global_iters = 2; % Number of MSEs to produce per shape.

%% SAVE MSE RESULTS TO FILE.
fid = fopen('Results_2d/mses.csv', 'wt');
% Set up for one global run.
% fprintf(fid, 'avg_mcmc_mse,avg_proj_mse,relative_change\n');
% Set up for multiple global runs on three shapes.
fprintf(fid, 'data_shape,gp,gp_proj,kern,kern_proj\n');


%% CONDUCT EXPERIMENT ON EACH SHAPE.
% List of shapes to run.
% shapes = {'trough', 'paraboloid', 'hand', 'parabolic_cylinder', ...
%           'wolverine', 'exponential', 'chair'};
% Try with only "flatter" surfaces.
shapes = {'hand', 'parabolic_cylinder', 'wolverine'};

% Run experiment for each shape.
for ii = 1:num_global_iters
    for shape = shapes
        shape_start_time = tic;
        
        [gp, gp_proj, kern, kern_proj] = gp_compare_mses_shape(tol_thres, ...
            eps1, eps2, iter, n, ls_factor, mesh_gran, num_posteriors, ...
            desired, d, shape{1}, fid);
        
        shape_time_elapsed = toc(shape_start_time);
        total_time_elapsed = toc(start_time);
        
        sendmail('momod@utexas.edu', 'UPDATE: MSE Comparisons', ...
            sprintf(strcat('Global iter: %s\n', ...
                           'Data shape:  %s\n', ...
                           'gp:                %s\n', ...
                           'gp_proj:        %s\n', ... 
                           'kern:             %s\n', ...
                           'kern_proj:     %s\n', ...
                           'shape_time:             %s\n', ...
                           'total_time_elapsed: %s\n'), ...
                num2str(ii), shape{1}, num2str(gp, '%0.2f'), ...
                num2str(gp_proj, '%0.2f'), num2str(kern, '%0.2f'), ...
                num2str(kern_proj, '%0.2f'), ...
                num2str(shape_time_elapsed, '%0.2f'), ...
                num2str(total_time_elapsed, '%0.2f')));
    end
end

fclose(fid);
   
sendmail('momod@utexas.edu', 'RESULTS: MSE Comparisons', '', ...
    {'/home/momod/Documents/gp/Results_2d/mses.csv'});

end_time = toc(start_time)



    
