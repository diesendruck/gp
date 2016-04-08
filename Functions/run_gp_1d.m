function [Eft_map] = run_gp_1d(x_nsy, y_nsy, ls_factor, x_grid)
% Run gpstuff, then produce and plot GP posterior curves.
%
% Inputs:
%   x_nsy: n x 1 matrix of data values.
%   y_nsy: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%   x_grid: Grid values over which to evaluate GP.
%
% Returns:
%   Eft_map: Sample from GP MAP.


%% STEP 0. Establish boundary of data, to make grid for curve.
x_range = max(x_nsy) - min(x_nsy);


%% STEP 1. Set up the GP.
noise_var_factor = 0.01;
length_scale = x_range*ls_factor;  % Scaled according to range.
mag_sig2 = (x_range*noise_var_factor)^2;  % Scaled according to range.

% Set up likelihood and covariance functions.
lik = lik_gaussian('sigma2', mag_sig2);
gpcf = gpcf_sexp('lengthScale', length_scale, 'magnSigma2', mag_sig2);

% Set up priors. Here, all parameters get uniform prior.
pn=prior_logunif();
lik = lik_gaussian(lik, 'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();

% Assemble covariance function with priors, and assemple gaussian process.
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);  % Assemble covariance function.
gp = gp_set('lik', lik, 'cf', gpcf);  % Assemble gaussian process.


%% STEP 2. Optimize GP and get params, MAP version.
opt = optimset('TolFun', 1e-3, 'TolX', 1e-3, 'Display','iter');
gp = gp_optim(gp, x_nsy, y_nsy, 'opt', opt);
[w,s] = gp_pak(gp);
disp(s), disp(exp(w))


%% STEP 3. Produce surface prediction.
[Eft_map, Varft_map] = gp_pred(gp, x_nsy, y_nsy, x_grid);


%% STEP 4. Plot the surface, with the training data.
% plot(x_grid, Eft_map, 'b-.');
% axis tight; hold on;
% plot(x_nsy, y_nsy, 'r.', 'Markers', 20);
% xlim([min(x_grid) max(x_grid)]);
% ylim([min(y_nsy)*1.1 max(y_nsy)*1.1]);

end

