function [xt1, xt2, z] = run_gp(x, y)
%   Run gp.

% GP REGRESSION ON TRAINING DATA ------------------------------------------
% Goal: Make GP regression data, finish with vectors x, y.

% Step 1. Train the GP.
lik = lik_gaussian('sigma2', 0.2^2);
gpcf = gpcf_sexp('lengthScale', [1.1 1.2], 'magnSigma2', 0.2^2)
pn=prior_logunif();
lik = lik_gaussian(lik, 'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gp = gp_set('lik', lik, 'cf', gpcf);

    %Should work now that X is multivariate
    %[K, C] = gp_trcov(gp, x);
    %example_x = [-1 -1 ; 0 0 ; 1 1];
    %[K, C] = gp_trcov(gp, example_x)

% Step 2. Optimize GP and get params, MAP version.
opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','iter');
gp=gp_optim(gp,x,y,'opt',opt);
[w,s] = gp_pak(gp);
disp(s), disp(exp(w))
% Step 2. Optimize GP and get params, MCMC version.  % --NOT WORKING YET--.
%[gp_rec,g,opt] = gp_mc(gp, x, y, 'nsamples', 220);
%gp_rec = thin(gp_rec, 21, 2);

% Step 3. Create surface grid.
[xt1,xt2]=meshgrid(-10.0:0.2:10.0,-10.0:0.2:10.0);
xt=[xt1(:) xt2(:)];

% Step 4. Produce surface prediction, MAP version.
[Eft_map, Varft_map] = gp_pred(gp, x, y, xt);
z = reshape(Eft_map, size(xt1));
% Step 4. Produce surface prediction, MCMC version.  % --NOT WORKING YET--.
%[Eft_s, Varft_s] = gpmc_preds(rfull, x, y, xt);
%[Eft_mc, Varft_mc] = gp_pred(gp_rec, x, y, xt);

% Step 4. Plot the surface, with the training data.
surf(xt1, xt2, z, 'FaceColor','interp', 'EdgeColor','flat', 'FaceLighting','gouraud');
axis tight; hold on;
plot3(x(:,1), x(:,2), y, 'r.', 'MarkerSize', 40);

end

