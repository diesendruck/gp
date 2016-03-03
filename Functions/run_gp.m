function [xt, Eft] = run_gp(x, y, ls_factor, mcmc_or_map)
% Run gpstuff to produce GP posterior surface.
% Inputs:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
%   ls_factor: Prior value for lengthscale hyperparameter.
%   color: Colormap selection, e.g. "winter" or "summer".

% STEP 0. Establish boundary of data, to make grid for surface.
grid_granularity = 20;
x_min = min(min(x)); x_max = max(max(x));  % Min/max across d dimensions.
x_range = x_max - x_min;
x_grid = x_min*1.1:x_range/grid_granularity:x_max*1.1;  % Create surface grid.
[xt1,xt2] = meshgrid(x_grid,x_grid);
xt=[xt1(:) xt2(:)];

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

% STEP 2. Optimize GP and get params, MAP version.
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
    [gp_rec, g, opt] = gp_mc(gp, x, y, 'nsamples', 120);
    gp_rec = thin(gp_rec, 21, 2);
    
    % STEP 3. Produce surface prediction.
    % [Eft_s, Varft_s] = gpmc_preds(rfull, x, y, xt);  % TODO: What is rfull?
    [Eft_mc, Varft_mc] = gp_pred(gp_rec, x, y, xt);
    z = reshape(Eft_mc, size(xt1));
    Eft = Eft_mc;

else
    error('Select either "MCMC" or "MAP" inference')

end

% STEP 4. Plot the surface, with the training data.
surf(xt1, xt2, z, 'FaceColor','interp', 'EdgeColor', 'flat', ...
    'FaceLighting','gouraud');
axis tight; hold on;
plot3(x(:,1), x(:,2), y, 'r.', 'MarkerSize', 40);

end

