% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html

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
n = 40;                % Data sample size.
d = 2                  % Dimension of data points.
ls_factor = 0.03;      % Lengthscale factor (proportion of x-range).
mesh_gran = 100;       % Number of ticks on mesh for plotting.
num_posteriors = 120;  % Number of posterior samples to generate.
desired = 3;           % Number of posterior samples to use.


%% CONDUCT EXPERIMENT ON EACH SHAPE.
% List of shapes to run.
shapes = {'trough', 'paraboloid', 'hand', 'parabolic_cylinder', ...
    'wolverine', 'exponential', 'chair'};

% Run experiment for each shape.
for shape = shapes
    gp_experiment_run_shape(tol_thres, eps1, eps2, iter, n, ls_factor, ...
        mesh_gran, num_posteriors, desired, d, shape{1});
end
    
    
    
    
    
    