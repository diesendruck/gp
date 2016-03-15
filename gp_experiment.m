% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html

tic

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
end

if 1  % Linux version.
    cd ~/Documents/gp/GPstuff-4.6/
    matlab_install
    cd ~/Documents/gp/Programs/
    run_mex_commands();
    cd ~/Documents/gp/
    addpath('~/Documents/gp/')
    addpath('~/Documents/gp/Programs')
    addpath('~/Documents/gp/Functions')
end

%% SET CONSTANTS.
tol_thres = 0;
eps1 = 10^-5;          % These 2 epsilons are used for convergence of the
eps2 = 10^-5;          %   convex projection algorithm.
iter = 0;              % Counter for iterations.
n = 30;                % Data sample size.
d = 2                  % Dimension of data points.
ls_factor = 0.03;      % Lengthscale factor (proportion of x-range).
mesh_gran = 30;       % Number of ticks on mesh for plotting.
num_posteriors = 1020;  % Number of posterior samples to generate.
desired = 100;           % Number of posterior samples to use.

%% SETUP EMAIL PARAMS.
myaddress = 'eltegog@gmail.com';
mypassword = 'T0g.eltegog';
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% SAVE MSE RESULTS TO FILE.
fid = fopen('mses.csv', 'wt');
fprintf(fid, 'avg_mcmc_mse,avg_proj_mse,relative_change\n');

%% CONDUCT EXPERIMENT ON EACH SHAPE.
% List of shapes to run.
shapes = {'trough', 'paraboloid', 'hand', 'parabolic_cylinder', ...
    'wolverine', 'exponential', 'chair'};

% Run experiment for each shape.
for shape = shapes
    gp_experiment_run_shape(tol_thres, eps1, eps2, iter, n, ls_factor, ...
        mesh_gran, num_posteriors, desired, d, shape{1}, fid);
    sendmail('momod@utexas.edu', 'UPDATES', sprintf('Finished %s', shape{1}));
end

fclose(fid);


sendmail('momod@utexas.edu', 'RESULTS', '', ...
    {'/home/momod/Documents/gp/mses.csv', ...
    '/home/momod/Documents/gp/trough.fig', ...
    '/home/momod/Documents/gp/trough_xq.csv', ...
    '/home/momod/Documents/gp/trough_yq.csv', ...
    '/home/momod/Documents/gp/trough_zq_proj.csv', ...  
    '/home/momod/Documents/gp/paraboloid.fig', ...
    '/home/momod/Documents/gp/paraboloid_xq.csv', ...
    '/home/momod/Documents/gp/paraboloid_yq.csv', ...
    '/home/momod/Documents/gp/paraboloid_zq_proj.csv', ...    
    '/home/momod/Documents/gp/hand.fig', ...
    '/home/momod/Documents/gp/hand_xq.csv', ...
    '/home/momod/Documents/gp/hand_yq.csv', ...
    '/home/momod/Documents/gp/hand_zq_proj.csv', ...
    '/home/momod/Documents/gp/parabolic_cylinder.fig', ...
    '/home/momod/Documents/gp/parabolic_cylinder_xq.csv', ...
    '/home/momod/Documents/gp/parabolic_cylinder_yq.csv', ...
    '/home/momod/Documents/gp/parabolic_cylinder_zq_proj.csv', ...
    '/home/momod/Documents/gp/wolverine.fig', ...
    '/home/momod/Documents/gp/wolverine_xq.csv', ...
    '/home/momod/Documents/gp/wolverine_yq.csv', ...
    '/home/momod/Documents/gp/wolverine_zq_proj.csv', ...
    '/home/momod/Documents/gp/exponential.fig', ...
    '/home/momod/Documents/gp/exponential_xq.csv', ...
    '/home/momod/Documents/gp/exponential_yq.csv', ...
    '/home/momod/Documents/gp/exponential_zq_proj.csv', ...
    '/home/momod/Documents/gp/chair.fig', ...
    '/home/momod/Documents/gp/chair_xq.csv', ...
    '/home/momod/Documents/gp/chair_yq.csv', ...
    '/home/momod/Documents/gp/chair_zq_proj.csv'});
    
toc



    
