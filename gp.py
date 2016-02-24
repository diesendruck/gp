import numpy as np
import GPy
from matplotlib import pyplot as pl
from __future__ import division

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
    model = GPy.models.GPRegression(X,Y,kernel, noise_var=1e-10)
    model.optimize(max_iters=10)
    #model.optimize('bfgs')

    testX = np.linspace(-10, 10, 100).reshape(-1, 1)
    posteriorTestY = model.posterior_samples_f(testX, full_cov=True, size=npost)
    simY, simMse = model.predict(testX)
    sigma = simMse**0.5

    fig = pl.figure()
    pl.plot(testX, simY, 'b-', label=u'Prediction')
    pl.fill(np.concatenate([testX, testX[::-1]]),
            np.concatenate([simY- 1.9600 * sigma,
                           (simY + 1.9600 * sigma)[::-1]]),
            alpha=.5, fc='b', ec='None', label='95% confidence interval')
    pl.errorbar(X.ravel(), Y.ravel(), dy, fmt='r.', markersize=10, label=u'Observations')
    pl.plot(testX, f(testX), 'r:', label=u'$f(x) = x^2$')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.legend(loc='upper center')
    pl.plot(testX, posteriorTestY)
    pl.plot(X, Y, 'ok', markersize=10)
    pl.show()

def main():
    n_data = 12
    n_post_samples = 10
    sample_GPs(n_data, n_post_samples)

main()

