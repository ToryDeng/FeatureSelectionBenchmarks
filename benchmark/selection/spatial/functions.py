# -*- coding: utf-8 -*-
# @Time : 2023/2/7 1:14
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
from multiprocessing import Process
from pathlib import Path
from typing import Union, Literal, Optional

import SpatialDE
import anndata as ad
import anndata2ri
import numpy as np
import pandas as pd
import scGeneClust as gc
from loguru import logger
from rpy2.robjects import pandas2ri, r, globalenv
from rpy2.robjects.packages import importr

from .._utils import HiddenPrints, select_features_from_df, rename_df_columns
from ..scrna.functions import SCRNASEQ_METHODS


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
    elif method == 'binspect_kmeans':
        selected_genes_df = BinSpect_kmeans(adata)
    elif method == 'binspect_rank':
        selected_genes_df = BinSpect_rank(adata)
    elif method in SCRNASEQ_METHODS.keys():
        selected_genes_df = SCRNASEQ_METHODS[method](adata)
    else:
        raise NotImplementedError(f"No implementation of {method}.")

    rename_df_columns(selected_genes_df)
    return select_features_from_df(selected_genes_df, n_selected_features)


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
        return pd.DataFrame({'Gene': raw_adata.var_names, 'Importance': 1 - results['qval']})


def SPARKX(adata: ad.AnnData) -> pd.DataFrame:
    with HiddenPrints():
        raw_adata = adata.raw.to_adata()
        pandas2ri.activate()
        spark = importr("SPARK")
        stats, res_stest, res_mtest = spark.sparkx(raw_adata.X.T, raw_adata.obsm['spatial'], verbose=False)
        pandas2ri.deactivate()
        results = pandas2ri.rpy2py(res_mtest).sort_values('adjustedPval')
        return pd.DataFrame({'Gene': raw_adata.var_names, 'Importance': 1 - results['adjustedPval']})


def Giotto_save_rds(adata: ad.AnnData):
    RDS_PATH = Path("./cache/temp")
    RDS_PATH.mkdir(parents=True, exist_ok=True)
    with HiddenPrints():
        anndata2ri.activate()
        importr("Giotto")
        globalenv['spe'] = adata
        r("""
        gobj = createGiottoObject(
        expression = assay(spe, 'X'),
        spatial_locs = reducedDim(spe, 'spatial'),
        cell_metadata = colData(spe),
        feat_metadata = rowData(spe),
        )
        saveRDS(gobj, './cache/temp/giotto_obj.rds')
        """)
        anndata2ri.deactivate()


def BinSpect_kmeans(adata: ad.AnnData) -> pd.DataFrame:
    raw_adata = adata.raw.to_adata()

    p = Process(target=Giotto_save_rds, args=(raw_adata,))
    p.start()
    p.join()
    with HiddenPrints():
        anndata2ri.activate()
        importr("Giotto")
        r("""
        gobj <- readRDS('./cache/temp/giotto_obj.rds')
    
        gobj <- normalizeGiotto(gobject = gobj, verbose = F)
        gobj <- createSpatialNetwork(gobject = gobj)
        """)
        result = r("binSpect(gobj, bin_method='kmeans', kmeans_algo = 'kmeans_arma_subset')")
        anndata2ri.deactivate()
    return pd.DataFrame({'Gene': result['feats'], 'Importance': 1 - result['adj.p.value']})


def BinSpect_rank(adata: ad.AnnData) -> pd.DataFrame:
    raw_adata = adata.raw.to_adata()

    p = Process(target=Giotto_save_rds, args=(raw_adata,))
    p.start()
    p.join()
    with HiddenPrints():
        anndata2ri.activate()
        importr("Giotto")
        r("""
        gobj <- readRDS('./cache/temp/giotto_obj.rds')
    
        gobj <- normalizeGiotto(gobject = gobj, verbose = F)
        gobj <- createSpatialNetwork(gobject = gobj)
        """)
        result = r("binSpect(gobj, bin_method='rank')")
        anndata2ri.deactivate()
    return pd.DataFrame({'Gene': result['feats'], 'Importance': 1 - result['adj.p.value']})
