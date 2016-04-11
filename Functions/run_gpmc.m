function [xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy,...
    y_nsy, ls_factor, num_posteriors, mesh_gran)
% Run gp mcmc and return samples of posterior.
%
% Args:
%   x_nsy: n x d matrix of data values.
%   y_nsy: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%   num_posteriors: Number of posterior samples to generate.
%   mesh_gran: Number of ticks on mesh for plotting.

%
% Returns:
%   xt: Matrix of grid points to evaluate over.
%   Eft_s: Samples from GP posterior.

do_diagnostics = 1;

%% STEP 0. Establish boundary of data, to make grid for surface.
[~, ~, ~, ~, x1_range, x2_range, xt1, xt2, xt] = compute_mesh_info(...
    x_nsy, mesh_gran);

%% STEP 1. Set up the GP.
noise_var_factor = 0.1;
length_scale = [x1_range*ls_factor, x2_range*ls_factor];  % Scaled according to range.
mag_sig2 = (min(x1_range, x2_range)*noise_var_factor)^2;  % Scaled according to range. % TODO: What should the sigma scaling be? MAX or MIN?

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
[rfull, g, opt] = gp_mc(gp, x_nsy, y_nsy, 'nsamples', num_posteriors);
gp_rec = thin(rfull, burned, thinned);
if do_diagnostics
    plot_mcmc_diagnostics(gp_rec, gp)
end

%% STEP 3. Produce surface prediction.
[Eft_s, Varft_s] = gpmc_preds(gp_rec, x_nsy, y_nsy, xt);  % Produce MCMC predictions.

% Print summary of results to console.
posterior_sample_count = size(Eft_s, 2);
sprintf('Of %d posterior samples, burned %d, thinned by %d: so %d remain.', ...
    [num_posteriors, burned, thinned, posterior_sample_count])

end

