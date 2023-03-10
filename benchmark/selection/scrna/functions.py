# -*- coding: utf-8 -*-
# @Time : 2023/2/6 1:31
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
import math
from typing import Optional, Union, Literal

import anndata as ad
import anndata2ri
import numpy as np
import pandas as pd
import scGeneClust as gc
import scanpy as sc
import triku as tk
from loguru import logger
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr

from ._giniclust3 import loessRegression, giniIndex
from .._utils import rpy2_wrapper, rename_df_columns, select_features_from_df


def select_in_single_batch(
        adata: ad.AnnData,
        method: str,
        n_selected_features: Union[int, Literal['auto']],
) -> np.ndarray:
    if method in SCRNASEQ_METHODS.keys():
        selected_genes_df = SCRNASEQ_METHODS[method](adata)
    else:
        raise NotImplementedError(f"No implementation of {method}.")

    rename_df_columns(selected_genes_df)
    # select features
    return select_features_from_df(selected_genes_df, n_selected_features)


def select_features(
        adata: ad.AnnData,
        method: str,
        n_selected_features: Union[int, Literal['auto']],
) -> np.ndarray:
    logger.opt(colors=True).info(
        f"Running <magenta>{method}</magenta> feature selection with "
        f"<yellow>{n_selected_features}</yellow> features to be selected."
    )
    if 'batch_key' in adata.uns_keys():
        batch_key = adata.uns['batch_key']
        if method.startswith('GeneClust'):
            raise ValueError(f"{method} is incompatible with the batch feature selection method in Seurat.")
        batch_genes = {}
        for batch in adata.obs[batch_key].unique():
            batch_adata = adata[adata.obs[batch_key] == batch, :]
            batch_genes[batch] = select_in_single_batch(batch_adata, method, n_selected_features)
        combined_genes = np.hstack(batch_genes.values())
        combined_ranks = np.hstack([np.arange(genes.shape[0]) for genes in batch_genes.values()])
        ugenes, ucounts = np.unique(combined_genes, return_counts=True)
        all_batch_df = pd.DataFrame({'n_datasets': ucounts}, index=ugenes)
        for gene in all_batch_df.index:
            all_batch_df.loc[gene, 'median_rank'] = np.median(combined_ranks[combined_genes == gene])
        all_batch_df.sort_values(by=['n_datasets', 'median_rank'], ascending=[False, True], inplace=True)
        selected_genes = all_batch_df.index.values[:n_selected_features]
    else:
        selected_genes = select_in_single_batch(adata, method, n_selected_features)
    return selected_genes


# python methods
def GeneClust_fast(adata: ad.AnnData) -> pd.DataFrame:
    selected_genes = gc.scGeneClust(adata.raw.to_adata(), n_var_clusters=200, version='fast', modality='sc')
    return pd.DataFrame({'Gene': selected_genes, 'Importance': np.ones(selected_genes.shape)})


def GeneClust_ps(adata: ad.AnnData) -> pd.DataFrame:
    n_obs_clusters = (~adata.obs[adata.uns['annot_key']].unique().isna()).sum().item()  # convert np.int64 to int
    logger.opt(colors=True).debug(f"GeneClust-ps will identify <yellow>{n_obs_clusters}</yellow> cell clusters.")
    selected_genes = gc.scGeneClust(adata.raw.to_adata(), n_obs_clusters=n_obs_clusters, version='ps', modality='sc')
    return pd.DataFrame({'Gene': selected_genes, 'Importance': np.ones(selected_genes.shape)})


def seurat_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=adata.n_vars // 2)
    return pd.DataFrame({'Gene': adata.var.dispersions_norm.index, 'Importance': adata.var.dispersions_norm})


def seurat_v3_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    raw_adata = adata.raw.to_adata()
    sc.pp.highly_variable_genes(raw_adata, flavor='seurat_v3', n_top_genes=adata.n_vars // 2)
    return pd.DataFrame({'Gene': raw_adata.var['variances_norm'].index, 'Importance': raw_adata.var['variances_norm']})


def triku_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    copied = adata.copy()  # .copy() prevent modification on adata object
    if 'X_pca' not in copied.obsm:
        sc.pp.pca(copied)
    if 'distances' not in copied.obsp or 'connectivities' not in copied.obsp:
        sc.pp.neighbors(copied)
    tk.tl.triku(copied, verbose='error')  # the larger the distance, the more important the gene is
    return pd.DataFrame({'Gene': copied.var_names, 'Importance': copied.var['triku_distance']})


def no_fs_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': np.ones(adata.var_names.shape)})


def giniclust3_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    geneGini, geneMax = giniIndex(adata.layers['normalized'].T)
    geneGini, geneLabel = np.array(geneGini), adata.var.index.tolist()
    logGeneMax = np.array(list(map(lambda x: math.log(x + 0.1, 2), geneMax)))
    _, sigGiniPvalue = loessRegression(geneGini, logGeneMax, geneLabel)
    result = pd.DataFrame.from_dict(sigGiniPvalue, orient='index', columns=['p.value']).reset_index()
    result['Importance'] = 1 - result['p.value']
    result.drop(columns=['p.value'], inplace=True)
    return result


def sc3_compute_importance(adata: ad.AnnData) -> pd.DataFrame:
    importance = np.zeros(shape=(adata.raw.shape[1],))  # the importances of filtered genes are set to zero
    dropouts = np.sum(adata.raw.X == 0, axis=0) / adata.raw.shape[1] * 100
    sc3_gene_filter = np.logical_and(10 < dropouts, dropouts < 90)
    logger.info(f"Filtering {adata.raw.shape[1] - sc3_gene_filter.sum()} rarely and ubiquitously expressed genes...")
    importance[sc3_gene_filter] = adata.raw.X[:, sc3_gene_filter].mean(axis=0)  # expression levels
    return pd.DataFrame({'Gene': adata.var_names, 'Importance': importance})


# R methods
@rpy2_wrapper
def FEAST(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    """
    Select features by FEAST. The input AnnData object contains norm and raw data. The raw data is normalized in R.

    Parameters
    ----------
    adata :
        anndata object.
    Returns
    -------
    A 1-column dataframe. The sole column contains genes sorted according to their F scores.
    The top genes are the most important.
    """
    importr('FEAST')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    Y <- process_Y(assay(sce, 'X'), thre = 2)
    Ynorm <- Norm_Y(Y)
    n_classes <- dim(unique(colData(sce)['celltype']))[1]
    rm(sce)
    idxs = FEAST_fast(Ynorm, k = n_classes)
    """)
    result = r("data.frame(Gene=rownames(Y)[idxs])")
    return result


@rpy2_wrapper
def M3Drop_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    """
    Select features by M3Drop. The input AnnData object contains norm and raw data. The raw data is normalized in R.

    Parameters
    ----------
    adata
      anndata object.
    Returns
    -------
    result
      A dataframe. The first column contains gene names, and the second column contains 1 - p.value.
    """
    importr('M3Drop')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    norm <- M3DropConvertData(assay(sce, 'X'), is.counts=TRUE)
    DE_genes <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=1, suppress.plot = TRUE)
    """)
    result = r("DE_genes").drop(columns=['effect.size', 'q.value'])
    result['Importance'] = 1 - result['p.value']
    result.drop(columns=['p.value'], inplace=True)
    return result


@rpy2_wrapper
def scmap_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('scmap')
    importr('scater')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    counts(sce) <- as.matrix(assay(sce, 'X'))
    sce <- logNormCounts(sce)
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- selectFeatures(sce, dim(sce)[1], suppress_plot = TRUE)
    res <- na.omit(rowData(sce)['scmap_scores'])
    """)
    result = r("res").reset_index()  # index  scmap_scores
    return result


@rpy2_wrapper
def deviance_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('scry')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())  # deviance need raw data
    result = r("rowData(devianceFeatureSelection(sce, assay='X'))['binomial_deviance']").reset_index()
    return result


@rpy2_wrapper
def scran(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('scuttle')
    importr('scran')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""
    assay(sce, 'X') <- as.matrix(assay(sce, 'X'))
    sce <- logNormCounts(sce, assay.type = "X")
    HVGs <- getTopHVGs(sce)
    """)
    result = r("data.frame(Gene=HVGs)")
    return result


@rpy2_wrapper
def sctransform_compute_importance(adata: ad.AnnData) -> Optional[pd.DataFrame]:
    importr('sctransform')
    globalenv['sce'] = anndata2ri.py2rpy(adata.raw.to_adata())
    r("""res <- vst(as(assay(sce, 'X'), "sparseMatrix"), method = "glmGamPoi", verbosity = 0)""")
    result = r("res$gene_attr['residual_variance']").reset_index()
    return result


SCRNASEQ_METHODS = {
    "GeneClust-fast": GeneClust_fast,
    'GeneClust-ps': GeneClust_ps,
    'seurat': seurat_compute_importance,
    'seurat_v3': seurat_v3_compute_importance,
    'triku': triku_compute_importance,
    'no_fs': no_fs_compute_importance,
    'giniclust3': giniclust3_compute_importance,
    'sc3': sc3_compute_importance,
    'feast': FEAST,
    'm3drop': M3Drop_compute_importance,
    'scmap': scmap_compute_importance,
    'deviance': deviance_compute_importance,
    'scran': scran,
    'sct': sctransform_compute_importance
}
