from __future__ import division
import numpy as np
import GPy
from matplotlib import pyplot as pl

def f(x):
    return x**2

def make_noisy_convex1(n):
    x = np.atleast_2d(np.random.sample(n)*20 - 10).T
    fx = f(x).ravel()
    dy = 0.5 + 5.0 * np.random.random(fx.shape)
    noise = np.random.normal(0, dy)
    y = fx + noise
    y = y.reshape(-1, 1)
    return x, y, dy

def sample_GPs(n_data, n_post_samples):
    n = n_data
    npost = n_post_samples
    X, Y, dy = make_noisy_convex1(n)

    kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
    m = GPy.models.GPRegression(X,Y,kernel, noise_var=1e-10)
    m.optimize(messages=True)
    m.optimize_restarts(num_restarts=10)

    testX = np.linspace(-10, 10, 100).reshape(-1, 1)
    posteriorTestY = m.posterior_samples_f(testX, full_cov=True, size=npost)
    simY, simMse = m.predict(testX)
    sigma = simMse**0.5

    fig = pl.figure(1)
    ax = fig.add_subplot(111)
    ax.fill(np.concatenate([testX, testX[::-1]]),
            np.concatenate([simY- 1.9600 * sigma,
                           (simY + 1.9600 * sigma)[::-1]]),
            alpha=.5, fc='b', ec='None', label='95% confidence interval')
    #pl.errorbar(X.ravel(), Y.ravel(), dy, fmt='r.', markersize=10, label=u'Observations')
    ax.plot(X.ravel(), Y.ravel(), 'r.', markersize=10, label=u'Observation')
    ax.plot(testX, f(testX), 'r--', label=u'$f(x) = x^2$')
    ax.plot(testX, simY, 'b-', label=u'Prediction')
    ax.plot(testX, posteriorTestY, ':')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1))
    fig.savefig('gp_posterior_sample', bbox_extra_artists=(lgd,),
                bbox_inches='tight', pad_inches=0.1)

    """
    fig = pl.figure(1)
    pl.fill(np.concatenate([testX, testX[::-1]]),
            np.concatenate([simY- 1.9600 * sigma,
                           (simY + 1.9600 * sigma)[::-1]]),
            alpha=.5, fc='b', ec='None', label='95% confidence interval')
    #pl.errorbar(X.ravel(), Y.ravel(), dy, fmt='r.', markersize=10, label=u'Observations')
    pl.plot(X.ravel(), Y.ravel(), 'r.', markersize=10, label=u'Observation')
    pl.plot(testX, f(testX), 'r--', label=u'$f(x) = x^2$')
    pl.plot(testX, simY, 'b-', label=u'Prediction')
    pl.plot(testX, posteriorTestY, ':', label=u'Posterior sample')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.legend(bbox_to_anchor=(0.5,0.75), loc='lower center', ncol=1)
    pl.show()
    """

def main():
    n_data = 5
    n_post_samples = 3
    sample_GPs(n_data, n_post_samples)

main()



    pl.show()
