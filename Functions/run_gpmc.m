function [xt1, xt2, xt, Eft_s, posterior_sample_count] = run_gpmc(x_nsy,...
    y_nsy, ls_factor, num_posteriors, mesh_gran)
% Run gp mcmc and return samples of posterior.
%
% Args:
%   x_nsy: n x d matrix of data values.
%   y_nsy: n x 1 matrix of response values.
%   ls_factor: Factor on prior value for lengthscale hyperparameter.
%   num_posteriors: Number of posterior samples to generate.
%   mesh_gran: Number of ticks on mesh for plotting.

%
% Returns:
%   xt: Matrix of grid points to evaluate over.
%   Eft_s: Samples from GP posterior.

% Plots MCMC traceplots. Turn off during global run to make figures work
% properly.
do_diagnostics = 0;


%% STEP 0. Establish boundary of data, to make grid for surface.
do_buffer = 0;
[~, ~, ~, ~, x1_range, x2_range, xt1, xt2, xt] = compute_mesh_info(...
    x_nsy, mesh_gran, do_buffer);


%% STEP 1. Set up the GP.
length_scale = [x1_range*ls_factor, x2_range*ls_factor];  % Scaled according to range.
mag_sig2 = 1;
lik_sig2 = 1;

% Set up likelihood and covariance functions.
gpcf = gpcf_sexp('lengthScale', length_scale, 'magnSigma2', mag_sig2);
lik = lik_gaussian('sigma2', lik_sig2);

% Set up priors.
%pmg = prior_sqrtunif();
%pns=prior_logunif();
%pls = prior_unif();
%pns = prior_sinvchi2('s2', 0.1,'nu', length(x_nsy));
%pns = prior_invgamma('sh', 4, 's', 1);
pmg = prior_invgamma('sh', 1, 's', 1);
pns = prior_sinvchi2('s2', 0.1,'nu', length(x_nsy));
pls = prior_invgamma('sh', 2, 's', length_scale(1));
pls = prior_gamma('sh', 15*length_scale(1), 'is', 5*length_scale(1));

% Assemble likelihood, covariance function, and gaussian process.
lik = lik_gaussian(lik, 'sigma2_prior', pns);
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pls, 'magnSigma2_prior', pmg);
gp = gp_set('lik', lik, 'cf', gpcf);


%% STEP 2. Optimize GP and get params.
burned = round(num_posteriors/5);
thinned = 3;
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

