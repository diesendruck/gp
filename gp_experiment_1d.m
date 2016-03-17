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

if 0  % Linux version.
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
n = 20;                % Data sample size.
d = 1;                  % Dimension of data points.
ls_factor = 0.01;      % Lengthscale factor (proportion of x-range).
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
fid = fopen('Results_1d/mses.csv', 'wt');
fprintf(fid, 'avg_mcmc_mse,avg_proj_mse,relative_change\n');

%% CONDUCT EXPERIMENT ON EACH SHAPE.
% List of shapes to run.
shapes = {'parabola'};

% Run experiment for each shape.
for shape = shapes
    gp_experiment_run_shape_1d(tol_thres, eps1, eps2, iter, n, ls_factor, ...
        mesh_gran, num_posteriors, desired, d, shape{1}, fid);
    sendmail('momod@utexas.edu', 'UPDATES', sprintf('Finished %s', shape{1}));
end

fclose(fid);

toc

sendmail('momod@utexas.edu', 'RESULTS_1D', '', ...
    {'/home/momod/Documents/gp/Results_1d/mses.csv', ...
    '/home/momod/Documents/gp/Results_1d/parabola.fig', ...
    '/home/momod/Documents/gp/Results_1d/parabola_x_grid.csv', ...
    '/home/momod/Documents/gp/Results_1d/parabola_avg_proj.csv', ...  
});
    



    
