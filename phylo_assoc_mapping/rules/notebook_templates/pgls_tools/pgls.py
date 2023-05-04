import numpy as  np
import pandas as pd

from statsmodels.regression.linear_model import GLS
from statsmodels.base.model import LikelihoodModelResults
from statsmodels.api import add_constant


def combine_results(results_list, model_type):
    # modified from statsmodels.imputation.mice.MICE.combine

    # Extract a few things from the models that were fit to
    # bootstrapped data sets.
    params_list = []
    cov_within = 0.
    scale_list = []
    for results in results_list:
        results_uw = results._results
        params_list.append(results_uw.params)
        cov_within += results_uw.cov_params()
        scale_list.append(results.scale)
    params_list = np.asarray(params_list)
    scale_list = np.asarray(scale_list)

    # The estimated parameters for the whole analysis
    params = params_list.mean(0)

    # The average of the within-bootstrap covariances
    cov_within /= len(results_list)

    # The between-bootstrap covariance
    cov_between = np.cov(params_list.T)

    # The estimated covariance matrix for the whole analysis
    f = 1 + 1 / float(len(results_list))
    cov_params = cov_within + f * cov_between

    # Set up a results instance
    scale = np.mean(scale_list)
    results = LikelihoodModelResults(
        model_type, params, cov_params / scale
    )
    results.scale = scale
    results.exog_names = results_list[0].model.exog_names
    results.endog_names = results_list[0].model.endog_names
    return results


def bootstrap_phylogenetic_gls(og, pheno, cov):
    bootstrap_res = []
    for _, pheno_boot in pheno.iterrows():
        pheno_boot = (pheno_boot - pheno_boot.mean()) / pheno_boot.std()
        model = GLS(endog=pheno_boot, exog=add_constant(og), sigma=cov)
        bootstrap_res.append(model.fit())
    res = combine_results(bootstrap_res, model_type=GLS)
    # get the 95% conf intervals for the coef
    coef = [r.params[1:] for r in bootstrap_res]
    lower, upper = np.percentile(coef, [2.5, 97.5], axis=0)
    return res.params[1:], lower, upper, res.pvalues[1:]