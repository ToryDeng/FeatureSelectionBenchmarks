# -*- coding: utf-8 -*-
# @Time : 2023/2/6 1:46
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import traceback

import anndata as ad
import anndata2ri
import numpy as np
import pandas as pd

from .._utils import HiddenPrints


def rpy2_wrapper(func):
    def wrapper(*args, **kwargs):
        try:
            with HiddenPrints():
                anndata2ri.activate()
                res = func(*args, **kwargs)
                anndata2ri.deactivate()
                return res
        except:
            traceback.print_exc()
    return wrapper


def select_features_from_df(selected_genes_df: pd.DataFrame, n_selected_features: int):
    if n_selected_features == 'auto':
        return selected_genes_df['Gene'].values
    else:
        if selected_genes_df.shape[0] < n_selected_features:
            msg = f"Genes are not enough to be selected ({selected_genes_df.shape[0]} exist(s), should select {n_selected_features})."
            raise RuntimeError(msg)
        else:
            return selected_genes_df.iloc[:n_selected_features, 0].values


def all_genes_in_adata(adata: ad.AnnData, selected_genes: np.ndarray):
    n_genes_in_adata = np.isin(selected_genes, adata.var_names).sum()
    if n_genes_in_adata != selected_genes.shape[0]:
        raise ValueError(f"Only {n_genes_in_adata} gene(s) are in adata, not {selected_genes.shape[0]}.")
    else:
        return True
