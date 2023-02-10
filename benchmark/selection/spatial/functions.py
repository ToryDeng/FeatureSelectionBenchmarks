# -*- coding: utf-8 -*-
# @Time : 2023/2/7 1:14
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
import os
from typing import Union, Literal, Optional

import SpatialDE
import anndata as ad
import numpy as np
import pandas as pd
import scGeneClust as gc
from loguru import logger
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

from .._utils import HiddenPrints, select_features_from_df, all_genes_in_adata


def select_features(
        adata: ad.AnnData,
        img: Optional[np.ndarray],
        method: str,
        n_selected_features: Union[int, Literal['auto']],
        random_state: int = 0
):
    logger.opt(colors=True).info(
        f"Running <magenta>{method}</magenta> feature selection with <yellow>{n_selected_features}</yellow> features to be selected.")
    if method == 'GeneClust-ps':
        selected_genes_df = GeneClust_ps(adata, img, random_state)
    elif method == 'spatialDE':
        selected_genes_df = spatialDE(adata)
    elif method == 'SPARKX':
        selected_genes_df = SPARKX(adata)
    else:
        raise NotImplementedError(f"No implementation of {method}.")

    selected_genes = select_features_from_df(selected_genes_df, n_selected_features)
    if all_genes_in_adata(adata, selected_genes):
        logger.info(f"Feature selection finished.")
        return selected_genes


def GeneClust_ps(adata: ad.AnnData, img: Optional[np.ndarray], random_state: int) -> pd.DataFrame:
    n_obs_clusters = (~adata.obs[adata.uns['annot_key']].unique().isna()).sum().item()  # convert np.int64 to int
    logger.opt(colors=True).debug(f"GeneClust-ps will identify <yellow>{n_obs_clusters}</yellow> spot clusters.")
    selected_genes = gc.scGeneClust(adata.raw.to_adata(), img, n_obs_clusters=n_obs_clusters, relevant_gene_pct=50, version='ps', modality='st', random_state=random_state)
    return pd.DataFrame({'Gene': selected_genes, 'Importance': np.ones(selected_genes.shape)})


def spatialDE(adata: ad.AnnData) -> pd.DataFrame:
    with HiddenPrints():
        raw_adata = adata.raw.to_adata()
        counts = raw_adata.to_df()
        coord = pd.DataFrame(raw_adata.obsm['spatial'], columns=['y_coord', 'x_coord'], index=raw_adata.obs_names)
        results = SpatialDE.run(coord, counts).sort_values('qval')
        return pd.DataFrame({'Gene': raw_adata.var_names, 'Importance': results['qval']})


def SPARKX(adata: ad.AnnData) -> pd.DataFrame:
    with HiddenPrints():
        raw_adata = adata.raw.to_adata()
        pandas2ri.activate()
        spark = importr("SPARK")
        stats, res_stest, res_mtest = spark.sparkx(raw_adata.X.T, raw_adata.obsm['spatial'], verbose=False)
        pandas2ri.deactivate()
        results = pandas2ri.rpy2py(res_mtest).sort_values('adjustedPval')
        return pd.DataFrame({'Gene': raw_adata.var_names, 'Importance': results['adjustedPval']})
