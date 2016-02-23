import numpy as np
from sklearn import gaussian_process
from matplotlib import pyplot as pl

# Create noisy version of convex, univariate data.
def make_noisy_convex():
    n_data_pts = 6
    x = np.atleast_2d(np.random.sample(n_data_pts)*20 - 10).T
    fx = f(x).ravel()
    dy = 0.5 + 3.0 * np.random.random(fx.shape)
    noise = np.random.normal(0, dy)
    y = fx + noise
    #pl.plot(x, y, 'r.')
    #pl.show()
    return x, y, dy

def f(x):
    return x**2

def fit_gp(X, Y, dy):
    gp = gaussian_process.GaussianProcess(verbose=True, corr='squared_exponential',
                                          theta0=1, thetaL=1e-3, thetaU=1)
    gp.fit(X, Y)

    n_mesh = 2000
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
            alpha=.25, fc='b', ec='None', label='95% confidence interval')
    pl.errorbar(X.ravel(), Y, dy, fmt='r.', markersize=10, label=u'Observations')
    pl.plot(mesh, f(mesh), 'r:', label=u'$f(x) = x^2$')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.ylim(-10, 110)
    pl.legend(loc='upper center')
    pl.show()
    return gp

def main():
    X, Y, dy = make_noisy_convex()
    posterior = fit_gp(X, Y, dy)
    posterior_cov = np.matrix(post.C) * np.matrix(post.C).T
    sample_posterior()


main()



