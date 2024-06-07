from typing import Sequence, Tuple, Collection, Union, Optional
import pandas as pd
import numpy as np
import datetime
import scipy
import scipy.linalg
import scipy.sparse
import os
import psutil
#
from scipy.sparse import csr_matrix
from anndata import AnnData
from scanpy import logging as logg
from anndata import logging as adlogg


def get_rss() -> str:
    return str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + 'GB'


def fmt_mem_usage_wdate(pfx_str) -> str:
    now = datetime.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    return current_time + ":" + pfx_str + ":" + adlogg.format_memory_usage(
        adlogg.get_memory_usage())


def log_mem_usage(phase) -> None:
    logg.info(fmt_mem_usage_wdate(phase))


def _aprior(delta_hat: np.ndarray) -> np.float_:
    m = delta_hat.mean()
    s2 = delta_hat.var()
    return (2 * s2 + m**2) / s2


def _bprior(delta_hat: np.ndarray) -> np.float_:
    m = delta_hat.mean()
    s2 = delta_hat.var()
    return (m * s2 + m**3) / s2


def sanitize_anndata(adata: AnnData) -> None:
    """Transform string annotations to categoricals."""
    adata._sanitize()


def _design_matrix(
    model: pd.DataFrame, batch_key: str, batch_levels: Collection[str]
) -> pd.DataFrame:
    """\
    Computes a simple design matrix.

    Parameters
    --------
    model
        Contains the batch annotation
    batch_key
        Name of the batch column
    batch_levels
        Levels of the batch annotation

    Returns
    --------
    The design matrix for the regression problem
    """
    import formulaic

    batch_formula: str = "~ 0 + C({}) ".format(batch_key)
    design: pd.DataFrame = formulaic.Formula(
        batch_formula).get_model_matrix(model, output="pandas")
    model: pd.DataFrame = model.drop([batch_key], axis=1)
    # TODO: numerical covariates not supported
    # numerical_covariates = model.select_dtypes('number').columns.values
    numerical_covariates: list[str] = []

    logg.info(f"-->Found {design.shape[1]} batches")
    other_cols: list[str] = [c for c in model.columns.values if c not in numerical_covariates]

    if other_cols:
        col_repr: str = " + ".join(x for x in other_cols)
        col_formula: str = "~ 0 + {}".format(col_repr)
        factor_matrix: pd.DataFrame = formulaic.Formula(
            col_formula).get_model_matrix(model, output="pandas")
        design: pd.DataFrame = pd.concat((design, factor_matrix), axis=1)
        logg.info(f"Found {len(other_cols)} categorical variables:")
        logg.info("\t" + ", ".join(other_cols) + '\n')

    log_mem_usage('Before leaving design matrix')
    return design


def row_sum_sparse(X_data: csr_matrix, y_data: np.ndarray,
                   ycfg_counts: np.ndarray,
                   ycfg_lookup: np.ndarray) -> np.ndarray:
    n_output: int = X_data.shape[0]
    in_dtype: np.dtype = X_data.dtype
    rsumdt: np.ndarray = np.zeros((n_output, 1), in_dtype)
    for bx, cfg_ctx in enumerate(ycfg_counts):
        cell_sel: np.ndarray = (ycfg_lookup == bx)
        ncells: np.ndarray = np.array(sum(cell_sel), dtype=in_dtype)
        if np.any(cell_sel):
            try:
                ydx: np.ndarray = (y_data[bx, :] * ncells).reshape((n_output, 1))
                ldx: np.ndarray = np.array(X_data[:, cell_sel].sum(axis=1),
                                           dtype=in_dtype)
                ldx -= ydx
                rsumdt += ldx.reshape((n_output, 1))
            except TypeError as ex:
                print(bx, cfg_ctx, sum(cell_sel), ex)
                raise
    return rsumdt


def row_mean_sparse(X_data: csr_matrix, y_data: np.ndarray,
                    ycfg_counts: np.ndarray,
                    ycfg_lookup: np.ndarray) -> np.ndarray:
    return np.divide(row_sum_sparse(X_data, y_data, ycfg_counts, ycfg_lookup),
                     np.array(X_data.shape[1], dtype=X_data.dtype))


def row_sum_of_squares_sparse(X_data: csr_matrix, y_data: np.ndarray,
                              ycfg_counts: np.ndarray,
                              ycfg_lookup: np.ndarray) -> np.ndarray:
    in_dtype: np.dtype = X_data.dtype
    ycfg_dcounts: np.ndarray = np.array(ycfg_counts, dtype=in_dtype)
    n_output: int = X_data.shape[0]
    xdsq: np.ndarray = np.array(X_data.multiply(X_data).sum(axis=1))
    # print(xdsq.shape)
    ydsq: np.ndarray = np.multiply(
        y_data ** 2,
        ycfg_dcounts.reshape(1, len(ycfg_dcounts))).sum(
            axis=1).reshape((n_output, 1))
    xydt: np.ndarray = np.zeros((n_output, 1), in_dtype)
    for bx, cfg_ctx in enumerate(ycfg_dcounts):
        cell_sel: np.ndarray = (ycfg_lookup == bx)
        if np.any(cell_sel):
            try:
                ydx: np.ndarray = y_data[:, bx].reshape(n_output, 1)
                xydt += np.array(X_data[:, cell_sel].multiply(ydx).sum(
                    axis=1), dtype=in_dtype)
            except TypeError as ex:
                print(bx, cfg_ctx, sum(cell_sel), ex)
                raise
    xydt = 2 * xydt
    # print(xdsq.shape, ydsq.shape, xydt.shape)
    return (xdsq + ydsq - xydt)


def _standardize_data(
    model: pd.DataFrame, data: csr_matrix, batch_key: str
) -> Tuple[csr_matrix, np.ndarray, pd.DataFrame, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray]:
    """\
    Standardizes the data per gene.

    The aim here is to make mean and variance be comparable across batches.

    Parameters
    --------
    model
        Contains the batch annotation
    data
        Contains the Data
    batch_key
        Name of the batch column in the model matrix

    Returns
    --------
    s_data
        Standardized Data
    design
        Batch assignment as one-hot encodings
    var_pooled
        Pooled variance per gene
    stand_mean
        Gene-wise mean
    """

    in_dtype: np.dtype = data.dtype
    # compute the design matrix
    batch_items = model.groupby(batch_key, observed=False).groups.items()
    batch_levels: Tuple[str]
    batch_info: Tuple[pd.Index]
    batch_levels, batch_info = zip(*batch_items)
    n_batch: int = len(batch_info)
    batch_counts: np.ndarray = np.array([len(v) for v in batch_info])
    n_array: float = float(sum(batch_counts))
    n_features: int = data.shape[0]

    design: pd.DataFrame = _design_matrix(model, batch_key, batch_levels)
    log_mem_usage('Obtained design matrix')
    # Assuming NO numerical covariates (TODO: numerical cov)
    # compute pooled variance estimator
    sp_design: csr_matrix = csr_matrix(design, dtype=in_dtype)
    cfg_vectors: np.ndarray
    cfg_lookup: np.ndarray
    cfg_counts: np.ndarray
    cfg_vectors, cfg_lookup, cfg_counts = np.unique(
        design, axis=0, return_counts=True, return_inverse=True)
    #
    B_hat: np.ndarray = np.dot(scipy.linalg.pinv((sp_design.T @ sp_design).A),
                               (sp_design.T @ data.T).A)
    # Mean only for the n_batch factors
    grand_mean: np.ndarray = np.dot(np.array(batch_counts / n_array, in_dtype).T,
                                    B_hat[:n_batch, :])
    # scanpy version:
    cfg_vectors: np.ndarray = np.array(cfg_vectors, dtype=in_dtype)
    cfgB_hat: np.ndarray = np.dot(cfg_vectors, B_hat)
    var_pooled: np.ndarray = row_sum_of_squares_sparse(data, cfgB_hat.T,
                                                       cfg_counts, cfg_lookup)
    var_pooled: np.ndarray = np.divide(var_pooled,
                                       n_array).reshape((n_features, 1))
    log_mem_usage('Done var_pooled')
    #
    # Compute the means
    if np.sum(var_pooled == 0) > 0:
        print(f'Found {np.sum(var_pooled == 0)} genes with zero variance.')
    tmp: np.ndarray = np.array(design.copy())
    tmp[:, :n_batch] = 0
    stdcfg_vecs: np.ndarray
    stdcfg_lookup: np.ndarray
    stdcfg_counts: np.ndarray
    stdcfg_vecs, stdcfg_lookup, stdcfg_counts = np.unique(tmp, axis=0,
                                                          return_counts=True,
                                                          return_inverse=True)
    stdcfg_vecs = np.array(stdcfg_vecs, dtype=in_dtype)
    standB_hat: np.ndarray = np.dot(stdcfg_vecs, B_hat)
    stand_mean: np.ndarray = standB_hat.T + grand_mean.reshape((n_features, 1))
    log_mem_usage('Done stand_mean')
    #
    # Right hand side - using broadcast
    s_dprime: csr_matrix = data.multiply(1 / np.sqrt(
        var_pooled.reshape((n_features, 1)))).tocsr()
    # Left hand side - element-wise division.
    s_mprime: np.ndarray = np.divide(
        stand_mean, np.sqrt(var_pooled.reshape((n_features, 1))))
    #
    # Need to be a bit careful with the zero variance genes
    # just set the zero variance genes to zero in the standardized data
    if np.any(var_pooled == 0):
        print("NVP", sum(var_pooled == 0))
        s_dprime[var_pooled == 0, :] = 0
        s_mprime[var_pooled == 0, :] = 0
    #
    log_mem_usage('About leave standardize')
    #
    Tuple[csr_matrix, np.ndarray, pd.DataFrame, np.ndarray, np.ndarray,
          np.ndarray, np.ndarray]
    return (s_dprime, s_mprime, design, var_pooled, stand_mean,  # type:ignore
            stdcfg_lookup, stdcfg_counts)


def _it_sol_sparse(s_dprime: csr_matrix, s_mprime: np.ndarray,
                   stdcfg_counts: np.ndarray, stdcfg_lookup: np.ndarray,
                   g_hat: np.ndarray, d_hat: np.ndarray,
                   g_bar: float, t2: float,
                   a: float, b: float, conv: float = 0.0001,
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """\
    Iteratively compute the conditional posterior means for gamma and delta.

    gamma is an estimator for the additive batch effect, deltat is an estimator
    for the multiplicative batch effect. We use an EB framework to estimate these
    two. Analytical expressions exist for both parameters, which however depend on each other.
    We therefore iteratively evalutate these two expressions until convergence is reached.

    Parameters
    --------
    s_data
        Contains the standardized Data
    g_hat
        Initial guess for gamma
    d_hat
        Initial guess for delta
    g_bar, t_2, a, b
        Hyperparameters
    conv: float, optional (default: `0.0001`)
        convergence criterium

    Returns:
    --------
    gamma
        estimated value for gamma
    delta
        estimated value for delta
    """

    # Assuming all NaN
    # n = (1 - np.isnan(s_dprime)).sum(axis=1)
    n = np.ones((s_dprime.shape[0],)) * s_dprime.shape[1]
    g_old = g_hat.copy()
    d_old = d_hat.copy()
    change = 1
    count = 0
    g_new = g_old
    d_new = d_old
    log_mem_usage('Before IT loop')
    # we place a normally distributed prior on gamma and and inverse gamma prior on delta
    # in the loop, gamma and delta are updated together. they depend on each other. we iterate until convergence.
    while change > conv:
        g_new = (t2 * n * g_hat + d_old * g_bar) / (t2 * n + d_old)
        # print(s_mprime.shape, g_new.shape)
        g_new_mp = s_mprime - g_new.reshape((g_new.shape[0], 1))
        sum2 = row_sum_of_squares_sparse(s_dprime, g_new_mp,
                                         stdcfg_counts, stdcfg_lookup)
        sum2 = sum2.flatten()
        d_new = (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)
        # print(sum2.shape, d_new.shape)
        #
        gsel = (g_old != 0)
        dsel = (d_old != 0)
        change = max((abs(g_new[gsel] - g_old[gsel]) / g_old[gsel]).max(),
                     (abs(d_new[dsel] - d_old[dsel]) / d_old[dsel]).max())
        g_old = g_new  # .copy()
        d_old = d_new  # .copy()
        count = count + 1
    log_mem_usage('After IT loop')

    return g_new, d_new


def sct_sparse(adata: AnnData, key: str = 'batch',
               covariates: Optional[Sequence[str]] = None,
               inplace: bool = True, run_dtype: type = np.float32,
               ) -> Union[AnnData, np.ndarray, None]:
    """\
    ComBat function for batch effect correction [Johnson07]_ [Leek12]_
    [Pedersen12]_.

    Corrects for batch effects by fitting linear models, gains statistical power
    via an EB framework where information is borrowed across genes.
    This uses the implementation `combat.py`_ [Pedersen12]_.

    .. _combat.py: https://github.com/brentp/combat.py

    Parameters
    ----------
    adata
        Annotated data matrix
    key
        Key to a categorical annotation from :attr:`~anndata.AnnData.obs`
        that will be used for batch effect removal.
    covariates
        Additional covariates besides the batch variable such as adjustment
        variables or biological condition. This parameter refers to the design
        matrix `X` in Equation 2.1 in [Johnson07]_ and to the `mod` argument in
        the original combat function in the sva R package.
        Note that not including covariates may introduce bias or lead to the
        removal of biological signal in unbalanced designs.
    inplace
        Whether to replace adata.X or to return the corrected data

    Returns
    -------
    Depending on the value of `inplace`, either returns the corrected matrix or
    or modifies `adata.X`.
    """

    # check the input
    if key not in adata.obs_keys():
        raise ValueError('Could not find the key {!r} in adata.obs'.format(key))

    if covariates is not None:
        cov_exist = np.isin(covariates, adata.obs_keys())
        if np.any(~cov_exist):
            missing_cov = np.array(covariates)[~cov_exist].tolist()
            raise ValueError(
                'Could not find the covariate(s) {!r} in adata.obs'.format(missing_cov)
            )

        if key in covariates:
            raise ValueError('Batch key and covariates cannot overlap')

        if len(covariates) != len(set(covariates)):
            raise ValueError('Covariates must be unique')
    else:
        covariates = []

    log_mem_usage('On Entry')
    # only works on dense matrices so far
    sanitize_anndata(adata)
    in_dtype = adata.X.dtype  # type: ignore
    if in_dtype == run_dtype:
        data: csr_matrix = adata.X.T  # type: ignore
    else:
        data: csr_matrix = adata.X.T.astype(run_dtype)  # type: ignore
    in_dtype = run_dtype
    log_mem_usage('After DF construct')

    # construct a pandas series of the batch annotation
    model: pd.DataFrame = adata.obs[[key] + list(covariates)]  # type:ignore
    batch_idxvals = model.groupby(key, observed=False).indices.values()
    n_batch = len(batch_idxvals)
    n_batches = np.array([len(v) for v in batch_idxvals])
    # n_array = float(sum(n_batches))
    n_features = data.shape[0]

    # standardize across genes using a pooled variance estimator
    logg.info("Standardizing Data across genes.\n")
    s_dprime, s_mprime, design, var_pooled, stand_mean, \
        stdcfg_lookup, stdcfg_counts = _standardize_data(model, data, key)
    log_mem_usage('Standardized Data')

    # fitting the parameters on the standardized data
    logg.info("Fitting L/S model and finding priors\n")
    batch_design: pd.DataFrame = design[design.columns[:n_batch]]  # type: ignore
    gamma_hat_pfx = np.linalg.pinv(np.diag(
        np.array(batch_design.T.sum(axis=1),
                 dtype=run_dtype))) @ batch_design.T
    sp_gamma_hat_pfx = scipy.sparse.csr_array(gamma_hat_pfx)
    # gamma_hat_dprime
    gamma_hat_dprime = (s_dprime @ sp_gamma_hat_pfx.T).A.T
    #
    n_stdcfgs = len(stdcfg_counts)
    batchxstd_counts = np.zeros((n_stdcfgs, n_batch), dtype=run_dtype)
    for jx in range(n_stdcfgs):
        batchxstd_counts[jx, :] = batch_design.loc[stdcfg_lookup == jx, :].sum(axis=0)
    # gamma_hat_mprime
    gamma_hat_mprime = (s_mprime @ batchxstd_counts) / n_batches
    gamma_hat = gamma_hat_dprime - gamma_hat_mprime.T
    log_mem_usage('Computed gamma_hat')

    # first estimate for the multiplicative batch effect
    delta_hat = []
    for i, batch_idxs in enumerate(batch_idxvals):
        # delta_b = np.zeros() TODO:
        cfg_bx = stdcfg_lookup[np.array(batch_design.iloc[:, i], dtype=bool)]
        batch_s_dprime = s_dprime[:, batch_idxs]
        mean_bx = s_mprime + row_mean_sparse(batch_s_dprime,
                                             s_mprime.T,
                                             batchxstd_counts[:, i],
                                             cfg_bx)
        var_bx = row_sum_of_squares_sparse(batch_s_dprime,
                                           mean_bx,
                                           batchxstd_counts[:, i],
                                           cfg_bx)
        #
        delta_hat.append(var_bx / (n_batches[i] - 1))
    log_mem_usage('Computed delta_hat')

    # empirically fix the prior hyperparameters
    gamma_bar = gamma_hat.mean(axis=1)
    t2 = gamma_hat.var(axis=1)
    # a_prior and b_prior are the priors on lambda and theta from Johnson and Li (2006)
    a_prior = list(map(_aprior, delta_hat))
    b_prior = list(map(_bprior, delta_hat))
    logg.info("--> Finding parametric adjustments")
    # -- TODO --
    # gamma star and delta star will be our empirical bayes (EB) estimators
    # for the additive and multiplicative batch effect per batch and cell
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_idxvals):
        # temp stores our estimates for the batch effect parameters.
        # temp[0] is the additive batch effect
        # temp[1] is the multiplicative batch effect
        gamma, delta = _it_sol_sparse(
            s_dprime[:, batch_idxs], s_mprime, batchxstd_counts[:, i],
            stdcfg_lookup[batch_idxs], gamma_hat[i], delta_hat[i].flatten(),
            gamma_bar[i], t2[i], float(a_prior[i]), float(b_prior[i]),)
        gamma_star.append(gamma)
        delta_star.append(delta)
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)
    log_mem_usage('ITR Complete All Batches')

    logg.info("--> Adjusting data")
    batch_adjustment = []
    # we now apply the parametric adjustment to the standardized data from above
    # loop over all batches in the data
    for j, batch_idxs in enumerate(batch_idxvals):
        # we basically substract the additive batch effect, rescale by the ratio
        # of multiplicative batch effect to pooled variance and add the overall gene
        # wise mean
        dsq = np.sqrt(delta_star[j, :])
        denom = dsq.reshape((len(dsq), 1))
        cfg_gamma_star = gamma_star[j, :].reshape((n_features, 1))
        batch_adjustment.append((cfg_gamma_star + s_mprime, denom))
    log_mem_usage('Compuated Batch Adjustments for all Batches')
    bayesdata = s_dprime.A
    log_mem_usage('Covert to data to Dense')
    for i in range(n_stdcfgs):
        for j, batch_idxs in enumerate(batch_idxvals):
            cell_sel = (stdcfg_lookup[batch_idxs] == i)
            numer_adj, denom = batch_adjustment[j]
            if np.any(cell_sel):
                bayesdata[:, batch_idxs[cell_sel]] -= numer_adj[:, i].reshape((
                    n_features, 1))
                bayesdata[:, batch_idxs[cell_sel]] /= denom

    log_mem_usage('Updated Batch Adjustments')
    vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
    # bayesdata = np.multiply(bayesdata, vpsq)
    bayesdata *= vpsq
    for i in range(len(stdcfg_counts)):
        bayesdata[:, stdcfg_lookup == i] += stand_mean[:, i].reshape(
            (n_features, 1))
    log_mem_usage('Normalized data')
    # put back into the adata object or return
    if inplace:
        adata.X = bayesdata.T
    else:
        return bayesdata.T
