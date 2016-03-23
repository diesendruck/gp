function [x_grid, Eft_s, posterior_sample_count] = run_gpmc_1d(x, ...
    y, ls_factor, num_posteriors)
% Run gp mcmc and return samples of posterior.
%
% Args:
%   x: n x 1 matrix of data values.
%   y: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%   num_posteriors: Number of posterior samples to generate.
%
% Returns:
%   x_grid: Array of grid points to evaluate over.
%   Eft_s: Samples from GP posterior.
%   posterior_sample_count: Posterior count after burn-in and thinning.

%% STEP 0. Establish boundary of data, to make grid for surface.
[x_l, x_h, x_range, x_grid] = compute_mesh_info_1d(x);

%% STEP 1. Set up the GP.
noise_var_factor = 0.01;
length_scale = x_range*ls_factor;  % Scaled according to range.
% Scaled according to range. 
% TODO: What should the sigma scaling be? MAX or MIN?
mag_sig2 = (min(x_range)*noise_var_factor)^2;  

% Set up likelihood and covariance functions.
lik = lik_gaussian('sigma2', mag_sig2);
gpcf = gpcf_sexp('lengthScale', length_scale, 'magnSigma2', mag_sig2);

% Set up priors. Here, all parameters get uniform prior.
pn=prior_logunif();
lik = lik_gaussian(lik, 'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();

% Assemble covariance function with priors, and assemple gaussian process.
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gp = gp_set('lik', lik, 'cf', gpcf);

%% STEP 2. Optimize GP and get params.
burned = 21;
thinned = 2;
[rfull, g, opt] = gp_mc(gp, x, y, 'nsamples', num_posteriors);
gp_rec = thin(rfull, burned, thinned);

%% STEP 3. Produce surface prediction.
[Eft_s, Varft_s] = gpmc_preds(gp_rec, x, y, x_grid);  % Produce MCMC predictions.

% Print summary of results to console.
posterior_sample_count = size(Eft_s, 2);
sprintf('Of %d posterior samples, burned %d, thinned by %d: so %d remain.', ...
    [num_posteriors, burned, thinned, posterior_sample_count])

end
