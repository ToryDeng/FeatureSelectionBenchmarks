# -*- coding: utf-8 -*-
# @Time : 2023/2/6 1:52
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
import anndata as ad
import anndata2ri
import numpy as np
import sc3s
import scanpy as sc
from loguru import logger
from rpy2.robjects import r, pandas2ri, globalenv
from rpy2.robjects.packages import importr
from sklearn.cluster import KMeans

from ..._utils import HiddenPrints


def cluster_cells(adata: ad.AnnData, method: str, k: int, random_state: int = 0):
    logger.opt(colors=True).info(f"Running <magenta>{method}</magenta> clustering with <yellow>{k}</yellow> clusters "
                                 f"and random state <yellow>{random_state}</yellow>...")
    if method == 'SHARP':
        return SHARP_clustering(adata, k, random_seed=random_state)
    elif method == 'Seurat_v4':
        logger.opt(colors=True).info(f"<magenta>{method}</magenta> clustering ignores the <yellow>`k={k}`</yellow> parameter.")
        return Seurat_v4_clustering(adata, random_seed=random_state)
    elif method == 'TSCAN':
        return TSCAN_clustering(adata, k)
    elif method == 'CIDR':
        return CIDR_clustering(adata, k)
    elif method == 'KMeans':
        return KMeans_clustering(adata, k, random_seed=random_state)
    elif method == 'SC3s':
        return SC3s_clustering(adata, k, random_seed=random_state)
    else:
        raise NotImplementedError(f"{method} has not been implemented.")


def SHARP_clustering(adata: ad.AnnData, k: int, random_seed: int):
    """
    Clustering cells using SHARP.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing raw counts and cell types
    k : int
        Number of clusters
    random_seed: int
        An integer for reproducibility


    Returns
    -------
    result : np.ndarray
        cluster labels
    """
    with HiddenPrints():
        anndata2ri.activate()
        importr('SHARP')
        globalenv['expr'], globalenv['k'] = pandas2ri.py2rpy(adata.to_df('log-normalized').T), k
        globalenv['seed'] = random_seed
        r("res <- SHARP(expr, N.cluster = k, prep = FALSE, rN.seed = seed, n.cores = 1)")
        label_pred = np.squeeze(np.array(r('res$pred_clusters')))
        anndata2ri.deactivate()
    return label_pred


def Seurat_v4_clustering(adata: ad.AnnData, random_seed: int):
    """
    Clustering cells using Seurat v4.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing raw counts and cell types
    random_seed: int
        An integer for reproducibility

    Returns
    -------
    result : np.ndarray
        cluster labels
    """
    assert adata.raw is not None, "The raw counts in data must exist!"
    raw_adata = adata.raw.to_adata()
    with HiddenPrints():
        anndata2ri.activate()
        importr('Seurat')
        importr('future')
        importr('doParallel')
        globalenv['sce'] = anndata2ri.py2rpy(raw_adata)
        globalenv['seed'] = random_seed
        r("""
        options(future.globals.maxSize= Inf)
        plan("multicore", workers = 6)
        as(sce, 'SingleCellExperiment')
        seuset <- as.Seurat(sce, counts='X', data=NULL)
        seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", verbose = FALSE)
        #seuset <- FindVariableFeatures(seuset)
        seuset <- ScaleData(object = seuset, verbose = FALSE)
        seuset <- RunPCA(object = seuset, features = rownames(seuset), verbose = FALSE)
        seuset <- FindNeighbors(object = seuset, verbose = FALSE)
        seuset <- FindClusters(object = seuset, random.seed = seed, verbose = FALSE)
        """)
        label_pred = np.squeeze(np.array(r('as.integer(unname(seuset$seurat_clusters))')))
        anndata2ri.deactivate()
    return label_pred


def SC3s_clustering(adata: ad.AnnData, k: int, random_seed: int):
    """
    Clustering cells using SC3s. Need normalized data.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing normalized data and cell types
    k : int
        Number of clusters
    random_seed: int
        An integer for reproducibility

    Returns
    -------
    result : np.ndarray
        cluster labels
    """
    with HiddenPrints():
        if 'X_pca' in adata.obsm:
            del adata.obsm['X_pca']  # do pca anyway
        sc.pp.pca(adata)
        cluster_adata = sc3s.tl.consensus(adata, n_clusters=[k], random_state=random_seed)
    pred_labels = cluster_adata.obs[f"sc3s_{k}"].to_numpy()
    del cluster_adata.obs[f"sc3s_{k}"]
    return pred_labels


def TSCAN_clustering(adata: ad.AnnData, k: int):
    with HiddenPrints():
        anndata2ri.activate()
        importr('TSCAN')
        globalenv['expr'], globalenv['k'] = pandas2ri.py2rpy(adata.to_df('log-normalized').T), k
        r("res <- exprmclust(expr, reduce = T, clusternum = k:k)")
        label_pred = np.squeeze(np.array(r('res$clusterid')))
        anndata2ri.deactivate()
    return label_pred


def KMeans_clustering(adata: ad.AnnData, k: int, random_seed: int):
    X_pca = sc.tl.pca(adata.X, return_info=False, random_state=random_seed)
    label_pred = KMeans(n_clusters=k, random_state=random_seed, n_init='auto').fit_predict(X_pca)
    return label_pred


def CIDR_clustering(adata: ad.AnnData, k: int):
    raw_adata = adata.raw.to_adata()
    with HiddenPrints():
        anndata2ri.activate()
        importr('cidr')
        globalenv['sce'], globalenv['k'] = anndata2ri.py2rpy(raw_adata), k
        r("""
        sData <- scDataConstructor(assay(sce, 'X'), tagType = "raw")
        sData <- determineDropoutCandidates(sData)
        sData <- wThreshold(sData)
        print('start to calculate dissimilarity...')
        sData <- scDissim(sData, threads = 0)
        sData <- scPCA(sData, plotPC = FALSE)
        sData <- nPC(sData)
        sDataC <- scCluster(object = sData, nCluster = k, nPC = sData@nPC, cMethod = "ward.D2")
        """)
        label_pred = np.squeeze(np.array(r("sDataC@clusters")))
        anndata2ri.deactivate()
    return label_pred
