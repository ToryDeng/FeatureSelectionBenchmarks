# -*- coding: utf-8 -*-
# @Time : 2023/2/7 19:17
# @Author : Tory Deng
# @File : _general.py
# @Software: PyCharm
from typing import Union, Callable, Literal, Optional

import anndata as ad
import numpy as np
import pandas as pd

from ._io import read_clusters_from_cache, write_clusters_as_cache
from .scrna.functions import cluster_cells
from .spatial.functions import cluster_spots


def generally_cluster_obs(
        fs_adata: ad.AnnData,
        img: Optional[np.ndarray],
        fs_method: Union[str, Callable],
        n_selected_genes: Union[int, Literal['auto']],
        cl_method: Union[str, Callable],
        run: int,
        modality: Literal['scrna', 'spatial'],
        **kwargs
):
    fs_method = fs_method.__name__ if callable(fs_method) else fs_method
    if isinstance(cl_method, str):
        # load cached clustering results or run clustering
        cluster_labels = read_clusters_from_cache(fs_adata.uns['data_name'], fs_method, n_selected_genes, cl_method, run)
        if cluster_labels is None:
            n_clusters = (~pd.isna(fs_adata.obs[fs_adata.uns['annot_key']].unique())).sum()
            if modality == 'scrna':
                cluster_labels = cluster_cells(fs_adata, cl_method, n_clusters, random_state=run)
            else:
                cluster_labels = cluster_spots(fs_adata, img, cl_method, n_clusters, random_state=run)
            write_clusters_as_cache(cluster_labels, fs_adata.uns['data_name'], fs_method, n_selected_genes, cl_method, run)
    elif callable(cl_method):  # cl_method is a function
        cluster_labels = read_clusters_from_cache(fs_adata.uns['data_name'], fs_method, n_selected_genes, cl_method.__name__, run)
        if cluster_labels is None:
            cluster_labels = cl_method(fs_adata, img, **kwargs)
            write_clusters_as_cache(cluster_labels, fs_adata.uns['data_name'], fs_method, n_selected_genes, cl_method.__name__, run)
    else:
        raise NotImplementedError(f"`cl_method` should be an valid string or a function, got {type(cl_method)}.")
    return cluster_labels

