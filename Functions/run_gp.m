function [xt, Eft] = run_gp(x, y, ls_factor, mcmc_or_map)
% Run gpstuff, then produce and plot GP posterior surface.
%
% Inputs:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%   mcmc_or_map: String "MCMC" or "MAP" to indicate which posterior to use.
%
% Returns:
%   xt: Matrix of grid points to evaluate over.
%   Eft: Sample from GP posterior.

%% STEP 0. Establish boundary of data, to make grid for surface.
[x1_l, x1_h, x2_l, x2_h, x1_range, x2_range, xt1, xt2, xt] = compute_mesh_info(x);

%% STEP 1. Set up the GP.
noise_var_factor = 0.01;
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
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);  % Assemble covariance function.
gp = gp_set('lik', lik, 'cf', gpcf);  % Assemble gaussian process.

%% STEP 2. Optimize GP and get params, MAP version.
if strcmp(mcmc_or_map, 'MAP')
    opt = optimset('TolFun', 1e-3, 'TolX', 1e-3, 'Display','iter');
    gp = gp_optim(gp, x, y, 'opt', opt);
    [w,s] = gp_pak(gp);
    disp(s), disp(exp(w))

    % STEP 3. Produce surface prediction.
    [Eft_map, Varft_map] = gp_pred(gp, x, y, xt);
    z = reshape(Eft_map, size(xt1));
    Eft = Eft_map;

% STEP 2. Optimize GP and get params, MCMC version.
elseif strcmp(mcmc_or_map, 'MCMC')
    [rfull, g, opt] = gp_mc(gp, x, y, 'nsamples', 120);
    gp_rec = thin(rfull, 21, 2);
    
    % STEP 3. Produce surface prediction.
    [Eft_s, Varft_s] = gpmc_preds(gp_rec, x, y, xt);
    % Sample one from GP MCMC posterior.
    num_samples = size(Eft_s, 2);
    Eft_smp = Eft_s(:, randi(num_samples));
    % Reshape output to enable graphing.
    z = reshape(Eft_smp, size(xt1));
    Eft = Eft_smp;

else
    error('Select either "MCMC" or "MAP" inference')

end

%% STEP 4. Plot the surface, with the training data.
surf(xt1, xt2, z, 'FaceColor','interp', 'EdgeColor', 'flat', ...
    'FaceLighting','gouraud');
axis tight; hold on;
plot3(x(:,1), x(:,2), y, 'r.', 'MarkerSize', 30);

end

