% This version was partially authored by Maurice Diesendruck.

% Adapted from Demo_example.m, a part of:
% Mazumder, R., Choudhury, A., Iyengar, G. and Sen, B. (2015).
%   A Computational Framework for Multivariate Convex Regression and its Variants.
%   Available at: http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html

% IMPORT GPSTUFF AND SET PATHS --------------------------------------------

cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/GPstuff-4.6/
matlab_install
addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/')
addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Programs')
cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/

% SET CONSTANTS -----------------------------------------------------------

tol_thres = 0;
eps1 = 10^-5;     %These 2 epsilons are used for convergence of the algo
eps2 = 10^-5;
iter = 0;         %counter for iterations
n = 10;            % Sample size
d = 2;            % Dimension d
sig = 7.0;        % Error variance

% SIMULATE RAW DATA (CONVEX + NOISE)
%[x, y] = make_noisy_convex(n, d, sig, 'trough');
[x, y] = make_noisy_convex(n, d, sig, 'paraboloid');
%scatter3(x(:, 1), x(:, 2), y);

% RUN GP ON RAW DATA.
figure;
subplot(1, 2, 1)
run_gp(x, y, 0.08, 'winter');  % Plots Gaussian Process surface.

% GET CONVEX PROJECTION
convex_y = project_to_convex(n, d, x, y, eps1, eps2, Max_Iter, rho);
%y-convex_y

% PLOT CONVEX OVER ORIGINAL GP
[xq, yq] = meshgrid(-10:.2:10);
vq = griddata(x(:,1), x(:,2), convex_y, xq, yq);
subplot(1,2,2)
mesh(xq,yq,vq);
hold on
plot3(x(:,1), x(:,2), convex_y, 'b.', 'MarkerSize', 40);

% Diagnostics from Sen, on outputs of project_to_convex that aren't
% currently being returned.
% h= figure(1)
% plot(cumsum(time_vec),sq_sse/n)
% title('Time versus Training SSE/n')     % SSE = sum of squared errors
% 
% h= figure(2)
% plot(cumsum(time_vec),prim_feas)
% title('Time versus Primal Feasibility/n')





