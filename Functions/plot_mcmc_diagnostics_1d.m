function plot_mcmc_diagnostics_1d(gp_rec, gp)
% Make traceplots and histogram of posterior hyperparameters.
%
% Args:
%   gp_rec: Thinned record of GP sampling from gpStuff.
%   gp: Full GP struct.
%
% Returns:
%   [none]

figure
clf

% Make traceplots and histogram for LS.
subplot(2, 3, 1)
plot(gp_rec.cf{1}.lengthScale(:,1));
title('Sample chain of length-scale')

subplot(2, 3, 4)
hist(gp_rec.cf{1}.lengthScale(:,1))
hold on
plot(gp.cf{1}.lengthScale(1), 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('Length-scale')

% Make traceplots and histogram for magnitude sigma2.
subplot(2, 3, 2)
plot(gp_rec.cf{1}.magnSigma2)
title('Sample chain of magnitude sig2')

subplot(2, 3, 5)
hist(gp_rec.cf{1}.magnSigma2)
hold on
plot(gp.cf{1}.magnSigma2, 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('Magnitude sig2')

% Make traceplots and histogram for noise variance.
subplot(2, 3, 3)
plot(gp_rec.lik.sigma2)
title('Sample chain of noise variance')

subplot(2, 3, 6)
hist(gp_rec.lik.sigma2)
hold on
plot(gp.lik.sigma2, 0, 'rx', 'MarkerSize', 11, 'LineWidth', 2)
title('Noise variance')


end

