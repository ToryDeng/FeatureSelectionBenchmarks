# -*- coding: utf-8 -*-
# @Time : 2023/2/8 18:51
# @Author : Tory Deng
# @File : test_benchmark.py
# @Software: PyCharm
from typing import Optional

import anndata as ad
import numpy as np

from benchmark._utils import rm_cache
from benchmark.run_benchmark import run_bench


def random_select(adata: ad.AnnData, n_selected_genes: int, seed: int):
    rng = np.random.default_rng(seed)
    return rng.choice(adata.var_names, size=n_selected_genes)


def random_clustering(adata: ad.AnnData, img: Optional[np.ndarray], k: int, seed: int):
    rng = np.random.default_rng(seed)
    return rng.integers(low=0, high=k, size=adata.n_obs)


def test_scrna_benchmark():
    data_cfg = {'PBMC3k': {
        'adata_path': 'tests/data/scrna/pbmc3k_raw.h5ad',
        'annot_key': 'louvain',
        'to_replace': {'Dendritic': ['Dendritic cells']},
        'to_remove': [np.nan, 'Megakaryocytes']
    }}
    fs_cfg = {'seurat_v3': [1000, 2000], random_select: [500], 'GeneClust-ps': ['auto']}
    cl_cfg = {'KMeans': 2}

    run_bench(data_cfg, fs_cfg, cl_cfg, ['ARI', 'NMI'], 'scrna', clean_cache=True, fs_kwarg={'seed': 123}, cl_kwarg={'k':5, 'seed': 123})
    rm_cache("./cache")


def test_spatial_benchmark():
    data_cfg = {'mouse_brain': {
        'adata_path': 'tests/data/spatial/V1_Adult_Mouse_Brain.h5ad',
        'image_path': 'tests/data/spatial/img.jpg',
        'annot_key': 'cluster',
    }}
    fs_cfg = {'spatialDE': [1000], 'SPARKX': [500], random_select: [2000], 'GeneClust-ps': ['auto']}
    cl_cfg = {'spaGCN': 1,  random_clustering: 1}

    run_bench(data_cfg, fs_cfg, cl_cfg, ['ARI', 'NMI'], 'spatial', clean_cache=True, fs_kwarg={'seed': 123}, cl_kwarg={'k':5, 'seed': 123})
    rm_cache("./cache")