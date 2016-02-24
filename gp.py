import numpy as np
import GPy
from sklearn import gaussian_process
import pyGPs
from matplotlib import pyplot as pl
from pylab import *
from __future__ import division
from numpy.linalg import inv

# Create noisy version of convex, univariate data.
def make_noisy_convex1():
    n_data_pts = 5
    x = np.atleast_2d(np.random.sample(n_data_pts)*20 - 10).T
    fx = f(x).ravel()
    dy = 0.5 + 5.0 * np.random.random(fx.shape)
    noise = np.random.normal(0, dy)
    y = fx + noise
    y = y[:,None]
    """
    n_data_pts = 15
    x = np.atleast_2d(np.random.sample(n_data_pts)*20 - 10).T
    dy = np.random.randn(*x.shape)*5.0 + 0.5
    y = x**2 + dy
    """
    return x, y, dy

def make_noisy_convex2():
    n_data_pts = 10
    x = np.atleast_2d(np.random.sample(n_data_pts)*10 - 5).T
    dy = np.random.randn(*x.shape)*1.0 + 0.5
    y = x**2 + dy
    return x, y, dy

def f(x):
    return x**2

def fit_gp1(X, Y, dy):
    Y = Y.ravel()
    gp = gaussian_process.GaussianProcess(verbose=True, corr='squared_exponential',
                                          theta0=1e-3, thetaL=1e-4, thetaU=1)
    gp.fit(X, Y)

    n_mesh = 100
    mesh = np.atleast_2d(np.linspace(-10, 10, n_mesh)).T
    y_pred, MSE = gp.predict(mesh, eval_MSE=True)
    sigma = np.sqrt(MSE)

    # Plot the function, the prediction and the 95% confidence interval based on
    # the MSE
    fig = pl.figure()
    pl.plot(mesh, y_pred, 'b-', label=u'Prediction')
    pl.fill(np.concatenate([mesh, mesh[::-1]]),
            np.concatenate([y_pred - 1.9600 * sigma,
                           (y_pred + 1.9600 * sigma)[::-1]]),
            alpha=.5, fc='b', ec='None', label='95% confidence interval')
    pl.errorbar(X.ravel(), Y.ravel(), dy, fmt='r.', markersize=10, label=u'Observations')
    pl.plot(mesh, f(mesh), 'r:', label=u'$f(x) = x^2$')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.legend(loc='upper center')
    pl.show()
    return gp, mesh, y_pred

def squared_exp(x1, x2, h):
    return np.exp(-0.5 * (x1-x2)**2 / h**2)

def fit_gp2(X, Y, dy):
    m = pyGPs.GPR()
    m.getPosterior(X, Y)
    m.setOptimizer("Minimize", num_restarts=10)
    m.optimize(X, Y)
    m.predict(mesh)
    m.plot()


    """
    # Make GP regression model with GPy.
    k1 = GPy.kern.RBF(1)
    m = GPy.models.GPRegression(X, Y)
    m.kern.lengthscale.set_prior(GPy.priors.Gamma.from_EV(1e-1, 5))
    m.kern.variance.set_prior(GPy.priors.Gamma.from_EV(1., 10.))
    m.constrain_positive('')
    m.optimize('bfgs')
    m.plot()
    print m

    n_mesh = 100
    mesh = np.atleast_2d(np.linspace(-10, 10, n_mesh)).T

    post_params = m.param_array
    post_k1 = GPy.kern.RBF(input_dim=1,
                       variance=post_params[0], lengthscale=post_params[1])
    post_k2 = GPy.kern.White(1, post_params[2])
    post_kernel = (post_k1 + post_k2).K(mesh)

    mu = m.predict(mesh)[0].ravel()
    plot(mesh, mu)
    plot(mesh, np.random.multivariate_normal(mu, post_kernel))
    plot(mesh, np.random.multivariate_normal(mu, m.kern.K(mesh)))

    plot(mesh, np.random.multivariate_normal(mu, np.atleast_2d(cscs)))

    mu, C = m.predict(mesh, full_cov=True)
    plot(mesh, np.random.multivariate_normal(mu.ravel(), C))


    hmc = GPy.inference.mcmc.HMC(m, stepsize=5e-2)
    s = hmc.sample(num_samples=1000)
    s = hmc.sample(num_samples=1000)
    plot(s)

    labels = ['kern variance', 'kern lengthscale','noise variance']
    samples = s[300:] # cut out the burn-in period
    from scipy import stats
    xmin = samples.min()
    xmax = samples.max()
    xs = np.linspace(xmin,xmax,100)
    for i in xrange(samples.shape[1]):
        kernel = stats.gaussian_kde(samples[:,i])
        plot(xs,kernel(xs),label=labels[i])

    _ = legend()
    return m, samples
    """

def sample_post(X, Y, mesh):
    h = 1
    noise_val = .3
    Kxx = squared_exp(X.ravel(), X, h)
    Kmx = squared_exp(mesh, X.ravel(), h)
    Kxm = Kmx.T
    Kmm = squared_exp(mesh.ravel(), mesh, h)
    noise_mat = np.eye(len(X))*noise_val
    post_mu = array(matrix(Kmx) * inv(Kxx + noise_mat) * Y).ravel()
    post_cov = Kmm - matrix(Kmx) * inv(Kxx + noise_mat) * Kxm

    plot(mesh, np.random.multivariate_normal(y_pred, Kmm))

def main():
    X, Y, dy = make_noisy_convex1()
    posterior1, mesh, y_pred = fit_gp1(X, Y, dy)
    chol = posterior1.C
    Kmm = matrix(chol) * matrix(chol).T
    Kmm = squared_exp(mesh.ravel(), mesh, h)
    plot(mesh, np.random.multivariate_normal(y_pred, Kmm))

main()
