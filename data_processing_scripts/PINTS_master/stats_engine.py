#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2019-08-07
#
# PINTS: Peak Identifier for Nascent Transcripts Sequencing
# Copyright (C) 2019 Li Yao at the Yu Lab
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
import sys
import warnings

try:
    import numpy as np
    from scipy.optimize import fmin_l_bfgs_b as BFGS
    from scipy.optimize import newton
    from scipy.special import digamma, factorial, gamma, gammaln, polygamma, psi
    from scipy.stats import binom_test, fisher_exact, nbinom, poisson

    # from statsmodels.base.model import GenericLikelihoodModel
except ImportError as e:
    missing_package = str(e).replace("No module named '", "").replace("'", "")
    sys.exit("Please install %s first!" % missing_package)
warnings.filterwarnings("error")


def zip_cdf(x, pi, lambda_):
    """

    Parameters
    ----------
    x
    pi
    lambda_

    Returns
    -------

    """
    assert x >= 0, "zip_cdf, x should > 0"
    p = pi + (1 - pi) * poisson.cdf(x, lambda_)
    if p > 1:
        p = 1
    elif p < 0:
        p = 0
    return p


def zip_moment_estimators(windows):
    """
    Moments estimators of ZIP
    :param windows: array-like
    :return: Corrected MMEs (lambda, pi)
    """
    s2 = windows.var()
    m = windows.mean()
    m2 = m**2
    if m >= s2:
        pi_mo = 0
        lamb_mo = m
    else:
        lamb_mo = (s2 + m2) / m - 1
        pi_mo = (s2 - m) / (s2 + m2 - m)

    return lamb_mo, pi_mo


def zip_probability_estimator(windows):
    """
    Estimate zip paras according to method mentioned in
    Zero-inflated models and estimation in zero-inflated Poisson distribution
    :param windows:
    :return: (lamb_mo, pi_pe)
    :ref: doi: 10.1080/03610918.2017.1341526
    """
    hat_p0 = sum(windows == 0) / windows.shape[0]
    s2 = windows.var()
    m = windows.mean()
    m2 = m**2
    if m >= s2:
        pi_pe = 0
        lamb_mo = m
    else:
        lamb_mo = (s2 + m2) / m - 1
        # pi_mo = (s2 - m) / (s2 + m2 - m)
        pi_pe = (hat_p0 - np.exp(-lamb_mo)) / (1 - np.exp(-lamb_mo))
    return lamb_mo, pi_pe


def zip_em(
    windows,
    init_lamb=None,
    init_pi=None,
    max_iter=1000,
    stop_diff=0.0001,
    debug=False,
    output_to="",
):
    """
    EM for Zero-inflated Poisson
    Parameters
    ----------
    windows : array-like

    init_lamb : float

    init_pi : float

    max_iter : int

    stop_diff : float

    debug : bool

    output_to : str

    Returns
    -------
    lambda_{mle} : float

    pi_{mle} : float

    is_convergent : bool

    llc : float
    """
    if init_lamb is None or init_pi is None:
        init_lamb, init_pi = zip_moment_estimators(windows=windows)
    lamb = init_lamb
    pi = init_pi
    n_iter = 0
    hat_z = np.zeros(len(windows))
    zero_elements = windows == 0

    if zero_elements.sum() == 0:
        return windows.mean(), 0, None, np.nan

    I = len(windows)
    prev_likelihood = 0.0
    likelihoods = []
    u_ele, c_ele = np.unique(windows, return_counts=True)
    while True:
        # expectation
        hat_z[zero_elements] = pi / (pi + np.exp(-lamb) * (1 - pi))
        # maximization
        pi = hat_z.sum() / I
        indicator = 1 - hat_z
        lamb = (indicator * windows).sum() / indicator.sum()
        # estimate likelihood
        u_pmf = (1 - pi) * poisson.pmf(u_ele, lamb)
        u_pmf[np.where(u_ele == 0)] += pi
        u_pmf[u_pmf < 10e-16] = 10e-16
        likelihood = (np.log(u_pmf) * c_ele).sum()
        likelihoods.append(likelihood)
        if n_iter > 0:
            if abs(likelihood - prev_likelihood) < stop_diff:
                if debug and output_to != "":
                    import matplotlib.pyplot as plt

                    plt.plot(likelihoods)
                    plt.ylabel("Log-likelihood")
                    plt.xlabel("Iteration")
                    plt.tight_layout()
                    plt.savefig(output_to, bbox_inches="tight", transparent=True)
                return lamb, pi, True, likelihood
            if n_iter > max_iter:
                return init_lamb, init_pi, False, likelihood
        prev_likelihood = likelihood
        n_iter += 1


def fit_nbinom(X, initial_params=None):
    infinitesimal = np.finfo(np.float).eps

    def log_likelihood(params, *args):
        r, p = params
        X = args[0]
        N = X.size

        # MLE estimate based on the formula on Wikipedia:
        # http://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
        result = (
            np.sum(gammaln(X + r))
            - np.sum(np.log(factorial(X)))
            - N * (gammaln(r))
            + N * r * np.log(p)
            + np.sum(X * np.log(1 - (p if p < 1 else 1 - infinitesimal)))
        )

        return -result

    if initial_params is None:
        # reasonable initial values (from fitdistr function in R)
        m = np.mean(X)
        v = np.var(X)
        size = (m**2) / (v - m) if v > m else 10

        # convert mu/size parameterization to prob/size
        p0 = size / ((size + m) if size + m != 0 else 1)
        r0 = size
        initial_params = np.array([r0, p0])

    try:
        bounds = [(infinitesimal, None), (infinitesimal, 1)]
        optimres = BFGS(
            log_likelihood,
            x0=initial_params,
            # fprime=log_likelihood_deriv,
            args=(X,),
            approx_grad=1,
            bounds=bounds,
        )

        params = optimres[0]
    except:
        # print("Failed to converge.")
        return r0, p0, False
    return params[0], params[1], True


def fit_nb(X, Z, pi, initial_params=None):
    infinitesimal = np.finfo(np.float).eps

    def log_likelihood(params, *args):
        mu, k = params
        X = args[0]
        Z = args[1]
        pi = args[2]
        muk = mu + k
        muk = muk if muk > 0 else infinitesimal
        pi = pi if pi > 0 else infinitesimal
        k = k if k > 0 else infinitesimal
        mu = mu if mu > 0 else infinitesimal
        one_minus_pi = 1 - pi
        one_minus_pi = one_minus_pi if one_minus_pi > 0 else infinitesimal
        # log_muk = np.log(muk if muk > 0 else infinitesimal)
        # N = len(X)
        is_zeros = X == 0
        is_zeros = is_zeros.astype(int)
        # incomplete log-likelihood function
        """
        result = np.sum(is_zeros * np.log(pi + (1 - pi) * ((k / (muk)) ** k))) + np.sum((1 - is_zeros) * (
                np.log(1 - pi) + gammaln(X + k) - gammaln(X + 1) - gammaln(k) + k * np.log(k) - k * np.log(
            muk) + X * np.log(mu if mu > 0 else infinitesimal) - X * np.log(muk)))
        """
        # complete log-likelihood function
        try:
            result = np.sum(Z * np.log(pi)) + np.sum(
                (1 - Z)
                * (
                    np.log(one_minus_pi)
                    + gammaln(X + k)
                    - gammaln(X + 1)
                    - gammaln(k)
                    + k * np.log(k)
                    - k * np.log(muk)
                    + X * np.log(mu)
                    - X * np.log(muk)
                )
            )
        except Exception:
            print("pi", pi, np.log(pi))
            print("1-pi", np.log(1 - pi))
            print(gammaln(k))
            print(np.log(k))
            print(np.log(mu))
            print(pi, mu, k, muk)
        return -result

    """
    def log_likelihood_deriv(params, *args):
        mu, k = params
        X = args[0]
        Z = args[1]
        one_minus_Z = 1 - Z
        mu_deriv = np.sum(one_minus_Z * X) / np.sum(one_minus_Z)
        k_deriv = np.sum(one_minus_Z * (psi(X + k) - psi(k) + np.log(k + 1) - np.log(mu + k) - ((k + X) / (mu + k))))
        return np.array([-mu_deriv, -k_deriv])
    """

    if initial_params is None:
        # reasonable initial values (from fitdistr function in R)
        m = np.mean(X)
        v = np.var(X)
        size = (m**2) / (v - m) if v > m else 10

        # convert mu/size parameters to mu/k
        p_0 = size / ((size + m) if size + m != 0 else 1)
        mu_0 = size * (1 - p_0) / p_0
        k_0 = 1 / size
        initial_params = np.array([mu_0, k_0])

    bounds = [(infinitesimal, None), (infinitesimal, None)]
    optimres = BFGS(
        log_likelihood,
        x0=initial_params,
        # fprime=log_likelihood_deriv,
        args=(X, Z, pi),
        approx_grad=1,
        bounds=bounds,
    )
    params = optimres[0]
    return params[0], params[1]


def zinb_em(
    windows,
    init_mu=None,
    init_k=None,
    init_pi=None,
    max_iter=1000,
    stop_diff=0.0001,
    debug=False,
    output_to="",
):
    """
    EM for Zero-inflated Negative Binomial
    :param windows:
    :param init_mu:
    :param init_k:
    :param init_pi:
    :param max_iter:
    :param stop_diff:
    :param debug:
    :param output_to:
    :return: (mu, k, pi, is_converged)
    """
    infinitesimal = np.finfo(np.float).eps
    if init_mu is None or init_k is None or init_pi is None:
        not_zero = windows != 0
        nn_zero = sum(not_zero)
        init_pi = (windows.shape[0] - not_zero.sum()) / windows.shape[0]
        if nn_zero > 0:
            init_mu = np.mean(windows[not_zero])
            s2 = np.var(windows[not_zero])
            size = init_mu**2 / (s2 - init_mu + 0.0001)
            size = size if size > 0 else 0.0001
            init_k = 1 / size
        else:
            init_mu = 0
            init_k = 1
    mu = init_mu
    k = init_k
    pi = init_pi
    mu_pre = mu
    k_pre = k
    pi_pre = pi
    n_iter = 0
    hat_z = np.zeros(len(windows))
    zero_elements = windows == 0
    n = len(windows)
    prev_likelihood = 0.0
    likelihoods = []

    while True:
        # expectation
        # in case of overflow
        # nb_term = k_pre / (mu_pre + k_pre)
        # if nb_term < 1 and k_pre >
        hat_z[zero_elements] = pi_pre / (
            pi_pre + (1 - pi_pre) * ((k_pre / (mu_pre + k_pre)) ** k_pre)
        )
        # maximization
        pi = hat_z.sum() / n
        # mu & k
        mu, k = fit_nb(windows, hat_z, pi_pre, [mu_pre, k_pre])
        # estimate likelihood
        pos_pmf = nbinom.pmf(windows, 1 / k, 1 / (1 + k * mu))
        pos_pmf[pos_pmf == 0] = infinitesimal
        likelihood = (
            np.sum(zero_elements * np.log(pi if pi > 0 else infinitesimal))
            + np.log(pos_pmf).sum()
        )
        likelihoods.append(likelihood)

        if n_iter > 0:
            if (
                abs(likelihood - prev_likelihood) < stop_diff
                or abs(mu - mu_pre) < stop_diff
                or abs(k - k_pre) < stop_diff
                or abs(pi - pi_pre) < stop_diff
            ):
                if debug and output_to != "":
                    import matplotlib.pyplot as plt

                    plt.plot(likelihoods)
                    plt.ylabel("Log-likelihood")
                    plt.xlabel("Iteration")
                    plt.tight_layout()
                    plt.savefig(output_to, bbox_inches="tight", transparent=True)
                return mu, k, pi, True, likelihood
            if n_iter > max_iter:
                return init_mu, init_k, init_pi, False, likelihood
        prev_likelihood = likelihood
        mu_pre = mu
        k_pre = k
        pi_pre = pi
        n_iter += 1


def zip_fit_mcmc(windows):
    """

    :param windows:
    :return:
    """
    import pymc3 as pm

    with pm.Model() as ZIP:
        psi = pm.Beta("p", 1, 1)
        lam = pm.Gamma("lam", 2, 0.1)

        y = pm.ZeroInflatedPoisson("y", lam, psi, observed=windows)
        trace = pm.sample(1000)
    pm.traceplot(trace[:])


def prop_test(pi_0, l_0, pi_1, l_1, empirical_threshold=5, alternative="greater"):
    """
    survivors = np.array([[1781, total1 - 1781], [1443, total2 - 47]])
    proportions_ztest
    In the two sample test, smaller means that the alternative hypothesis is
    p1 < p2 and larger means p1 > p2 where p1 is the proportion of the first sample and p2 of the second one
    :return:
    """
    count_0 = int(pi_0 * l_0)
    total_0 = l_0
    count_1 = int(pi_1 * l_1)
    total_1 = l_1
    # if no zero
    if count_0 == 0 or count_1 == 0:
        return 10e-16
    try:
        while count_0 < empirical_threshold:
            count_0 = int(pi_0 * l_0 * 10)
            total_0 *= 10
        while count_1 < empirical_threshold:
            count_1 = int(pi_1 * l_1 * 10)
            total_1 *= 10
        _, pval = fisher_exact(
            np.array([[count_0, count_1], [total_0 - count_0, total_1 - count_1]]),
            alternative=alternative,
        )
    except Exception as e:
        print(e, count_0, total_0, count_1, total_1)
        pval = 1
    return pval


def poisson_test(observed_lambda, background_lambda):
    observed_scale = 1
    background_scale = 1
    while observed_lambda < 1:
        observed_scale *= 10
        observed_lambda *= observed_scale
    while background_lambda < 1:
        background_scale *= 10
        background_lambda *= background_scale
    observed_total = int(10000 * observed_scale)
    background_total = int(background_scale * 10000)
    return binom_test(
        (int(observed_lambda), int(background_lambda)),
        observed_total / (observed_total + background_total),
        alternative="less",
    )


def get_outlier_threshold(data, direction=-1):
    """
    Get outlier threshold

    Parameters
    ----------
    data :
    direction : int
        Which direction of outlier threshold to be returned. 1 for upper, 0 for both, -1 for lower
    Returns
    -------
    lower_bound : int or None
        Lower threshold for outlier detection
    upper_bound : int or None
        Upper threshold for outlier detection
    """
    q1 = np.quantile(data, 0.25)
    q3 = np.quantile(data, 0.75)
    iqr = q3 - q1
    if direction == 1:
        return None, q3 + 1.5 * iqr
    elif direction == -1:
        return q1 - 1.5 * iqr, None
    else:
        return q1 - 1.5 * iqr, q3 + 1.5 * iqr


def pval_dist(pval_list, logger, output_diagnostic_plot=True, output_to=None):
    """
    Plot pval distribution (histogram)

    Parameters
    ----------
    pval_list : array-like
        p-values
    logger : Python Logger Objects
        Logger for recording errors / info
    output_diagnostic_plot : bool
        Whether to generate the histogram
    output_to : str or None
        Path for the output file

    Returns
    -------

    """
    pval_list_filtered = pval_list[np.logical_and(pval_list >= 0.0, pval_list <= 1)]
    if output_diagnostic_plot and output_to is not None:
        if len(pval_list_filtered) < 20:
            logger.error("Cannot get enough p-values for binning")
            return
        import matplotlib.pyplot as plt

        plt.hist(
            pval_list_filtered,
            bins=20,
            range=(0, 1),
            density=True,
            color="powderblue",
        )
        plt.xlabel("$p$-value")
        plt.ylabel("Density")
        plt.xlim((0, 1))
        plt.tight_layout()
        plt.savefig(output_to, transparent=True, bbox_inches="tight")
        plt.close()
        logger.info("Diagnostic plot for p-values was wrote to %s" % output_to)


if __name__ == "__main__":
    np.random.seed(299)
    n = 1000
    theta = 2.5  # Poisson rate
    pi = 0.55  # probability of extra-zeros (pi = 1-psi)
    mu = 4.48
    k = 0.25

    # Simulate some data
    counts = np.array(
        [(np.random.random() > pi) * np.random.poisson(theta) for i in range(n)]
    )
    print(zip_em(counts, debug=True, output_to="zip_ll.pdf"))
