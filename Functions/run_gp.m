function [xt1, xt2, z] = run_gp(x, y, ls_factor, color)
% Run gpstuff.

% GP REGRESSION ON TRAINING DATA ------------------------------------------
% Goal: Make GP regression data, finish with vectors x, y.

% Establish boundary of data.
x_min = min(min(x)); x_max = max(max(x)); y_min = min(y); y_max = max(y);
x_range = x_max - x_min;

% Step 1. Train the GP.
length_scale = [x_range*ls_factor, x_range*ls_factor];  % Scaled according to range.
mag_sig2 = (x_range*0.01)^2;

lik = lik_gaussian('sigma2', mag_sig2);
gpcf = gpcf_sexp('lengthScale', length_scale, 'magnSigma2', mag_sig2);

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
x_grid = x_min*1.1:x_range/100:x_max*1.1;
[xt1,xt2]=meshgrid(x_grid,x_grid);  %Scaled according to range.
xt=[xt1(:) xt2(:)];

% Step 4. Produce surface prediction, MAP version.
[Eft_map, Varft_map] = gp_pred(gp, x, y, xt);
z = reshape(Eft_map, size(xt1));
% Step 4. Produce surface prediction, MCMC version.  % --NOT WORKING YET--.
%[Eft_s, Varft_s] = gpmc_preds(rfull, x, y, xt);
%[Eft_mc, Varft_mc] = gp_pred(gp_rec, x, y, xt);

% Step 4. Plot the surface, with the training data.
surf(xt1, xt2, z, 'FaceColor','interp', 'EdgeColor','flat', 'FaceLighting','gouraud');
if color=='winter'
    colormap winter;
else
    colormap summer;
end
axis tight; hold on;
plot3(x(:,1), x(:,2), y, 'r.', 'MarkerSize', 40);

end

