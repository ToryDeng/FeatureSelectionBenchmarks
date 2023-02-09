# -*- coding: utf-8 -*-
# @Time : 2023/2/7 10:37
# @Author : Tory Deng
# @File : functions.py
# @Software: PyCharm
import contextlib
import os
import random

import SpaGCN as spg
import anndata as ad
import numpy as np
import torch

from ..._utils import HiddenPrints

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
import stlearn as st
from pathlib import Path
from loguru import logger


def cluster_spots(adata: ad.AnnData, image: np.ndarray, method: str, k: int, random_state: int = 0):
    logger.opt(colors=True).info(f"Running <magenta>{method}</magenta> clustering with <yellow>{k}</yellow> clusters "
                                 f"and random state <yellow>{random_state}</yellow>...")
    if method == 'spaGCN':
        return spaGCN_clustering(adata, image, k, random_state=random_state)
    elif method == 'stlearn':
        return stlearn_clustering(adata, k, random_state=random_state)
    else:
        raise NotImplementedError(f"{method} has not been implemented.")


def spaGCN_clustering(
        adata: ad.AnnData,
        img: np.ndarray,
        k: int,
        # shape: Literal['hexagon', 'square'] = 'hexagon',
        random_state: int = 0
):
    # prepare positional information
    x_array, y_array = adata.obs["array_row"].values, adata.obs["array_col"].values
    x_pixel, y_pixel = adata.obsm['spatial'][:, 1], adata.obsm['spatial'][:, 0]

    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
        if img is None:
            adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)
        else:
            adj = spg.calculate_adj_matrix(
                x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=49, alpha=1, histology=True
            )
        l = spg.search_l(0.5, adj, start=0.01, end=1000, tol=0.01, max_run=100)
        res = spg.search_res(adata, adj, l, k, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
                             r_seed=random_state, t_seed=random_state, n_seed=random_state)
        clf = spg.SpaGCN()
        clf.set_l(l)
        random.seed(random_state)
        torch.manual_seed(random_state)
        np.random.seed(random_state)
        clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob = clf.predict()
        shape = adata.uns['shape'] if 'shape' in adata.uns else 'hexagon'
        return np.array(spg.spatial_domains_refinement_ez_mode(adata.obs.index, y_pred, x_array, y_array, shape))


def stlearn_clustering(adata: ad.AnnData, k: int, random_state: int = 0):
    with HiddenPrints():
        TILE_PATH = Path("./cache/temp/stlearn/tiles")
        TILE_PATH.mkdir(parents=True, exist_ok=True)

        copied_adata = adata.copy()
        library_id = list(copied_adata.uns["spatial"].keys())[0]
        copied_adata.uns["spatial"][library_id]["use_quality"] = 'hires'
        copied_adata.obs['imagerow'], copied_adata.obs['imagecol'] = copied_adata.obs['array_row'], copied_adata.obs['array_col']

        st.em.run_pca(copied_adata, random_state=random_state, use_highly_variable=False)
        st.pp.tiling(copied_adata, TILE_PATH)
        st.pp.extract_feature(copied_adata, seeds=random_state)
        st.spatial.SME.SME_normalize(copied_adata, use_data='raw', weights="physical_distance")
        copied_adata.X = copied_adata.obsm['raw_SME_normalized']
        st.pp.scale(copied_adata)
        st.em.run_pca(copied_adata, random_state=random_state)
        st.tl.clustering.kmeans(copied_adata, n_clusters=k, use_data="X_pca", key_added="X_pca_kmeans")

        return copied_adata.obs['X_pca_kmeans'].values
