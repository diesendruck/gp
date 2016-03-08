function [xt, Eft_s] = run_gpmc(x, y, ls_factor)
% Run gp mcmc and return samples of posterior.
%
% Inputs:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%
% Returns:
%   xt: Matrix of grid points to evaluate over.
%   Eft_s: Samples from GP posterior.

% STEP 0. Establish boundary of data, to make grid for surface.
[xt, x_range] = compute_meshgrid_matrix(x);

% STEP 1. Train the GP.
noise_var_factor = 0.01;
length_scale = [x_range*ls_factor, x_range*ls_factor];  % Scaled according to range.
mag_sig2 = (x_range*noise_var_factor)^2;  % Scaled according to range.

lik = lik_gaussian('sigma2', mag_sig2);
gpcf = gpcf_sexp('lengthScale', length_scale, 'magnSigma2', mag_sig2);

pn=prior_logunif();  % All parameters get uniform prior.
lik = lik_gaussian(lik, 'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();

gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);  % Assemble covariance function.
gp = gp_set('lik', lik, 'cf', gpcf);  % Assemble gaussian process.

% STEP 2. Optimize GP and get params.
[rfull, g, opt] = gp_mc(gp, x, y, 'nsamples', 120);
gp_rec = thin(rfull, 21, 2);
    
% STEP 3. Produce surface prediction.
[Eft_s, Varft_s] = gpmc_preds(gp_rec, x, y, xt);  % Produce MCMC predictions.

end

