function [Eft_s, posterior_sample_count] = run_gpmc_1d(x_nsy, y_nsy, ...
    ls_factor, num_posteriors, mesh_gran)
% Run gp mcmc and return samples of posterior.
%
% Args:
%   x_nsy: n x 1 matrix of data values.
%   y_nsy: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%   num_posteriors: Number of posterior samples to generate.
%   mesh_gran: Number of ticks on mesh for plotting.
%
% Returns:
%   Eft_s: Samples from GP posterior.
%   posterior_sample_count: Posterior count after burn-in and thinning.

% Plots MCMC traceplots. Turn off during global run to make figures work
% properly.
do_diagnostics = 0;


%% STEP 0. Establish boundary of data, to make grid for surface.
[x_l, x_h, x_range, x_grid] = compute_mesh_info_1d(x_nsy, mesh_gran);


%% STEP 1. Set up the GP.
length_scale = x_range*ls_factor;  % Scaled according to range.
mag_sig2 = 1;
lik_sig2 = 1;

% Set up likelihood and covariance functions.
gpcf = gpcf_sexp('lengthScale', length_scale, 'magnSigma2', mag_sig2);
lik = lik_gaussian('sigma2', lik_sig2);

% Set up priors.
%pn=prior_logunif();
%pl = prior_unif();
%pm = prior_sqrtunif();
%pmg = prior_invgamma('sh', 1, 's', 1);
pmg = prior_gamma('sh', 10, 'is', 10);
pns = prior_sinvchi2('s2', 0.1,'nu', length(x_nsy));
%pls = prior_invgamma('sh', 2, 's', length_scale(1));
pls = prior_gamma('sh', 15*length_scale, 'is', 5*length_scale);

% Assemble covariance function with priors, and assemple gaussian process.
lik = lik_gaussian(lik, 'sigma2_prior', pns);
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pls, 'magnSigma2_prior', pmg);
gp = gp_set('lik', lik, 'cf', gpcf);


%% STEP 2. Optimize GP and get params.
burned = round(num_posteriors/4);
thinned = 4;
[rfull, g, opt] = gp_mc(gp, x_nsy, y_nsy, 'nsamples', num_posteriors);
gp_rec = thin(rfull, burned, thinned);
if do_diagnostics
    plot_mcmc_diagnostics_1d(gp_rec, gp)
end


%% STEP 3. Produce surface prediction.
[Eft_s, Varft_s] = gpmc_preds(gp_rec, x_nsy, y_nsy, x_grid);  % Produce MCMC predictions.

% Print summary of results to console.
posterior_sample_count = size(Eft_s, 2);
sprintf('Of %d posterior samples, burned %d, thinned by %d: so %d remain.', ...
    [num_posteriors, burned, thinned, posterior_sample_count])


end

