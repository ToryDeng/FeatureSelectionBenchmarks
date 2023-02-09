# -*- coding: utf-8 -*-
# @Time : 2023/2/6 1:06
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import anndata as ad
import numpy as np
import scanpy as sc
from loguru import logger
from scipy.sparse import issparse


def subsample(adata: ad.AnnData, n_sampled_obs: int = None, n_sampled_vars: int = None, random_state: int = 0):
    assert n_sampled_obs is not None or n_sampled_vars is not None, \
        "Specify at least one of `n_sampled_obs` and `n_sampled_vars`."
    if n_sampled_obs is not None:
        adata = sc.pp.subsample(adata, n_obs=n_sampled_obs, random_state=random_state, copy=True)
        sc.pp.filter_genes(adata, min_cells=1)  # filter genes with no count

    if n_sampled_vars is not None:
        rng = np.random.default_rng(random_state)
        adata = adata[:, rng.choice(adata.n_vars, size=int(n_sampled_vars), replace=False)].copy()
        sc.pp.filter_cells(adata, min_genes=1)  # filter cells with no count

    logger.debug(f"Shape of sampled adata : {adata.shape}")


def to_dense(adata: ad.AnnData):
    if issparse(adata.X):
        logger.info(f"Found sparse matrix in `adata.X`. Converting to ndarray...")
        adata.X = adata.X.toarray()


def is_normalized(adata: ad.AnnData):
    return not np.allclose(adata.X % 1, 0)
