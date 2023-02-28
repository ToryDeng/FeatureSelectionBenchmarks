# -*- coding: utf-8 -*-
# @Time : 2023/2/6 2:08
# @Author : Tory Deng
# @File : main.py
# @Software: PyCharm


# TODO: add Seurat clustering for SRT
# TODO: check GeneClust dependency

from benchmark._utils import rm_cache
from benchmark.run_benchmark import run_bench

data_cfg = {'mouse_brain': {
        'adata_path': 'tests/data/spatial/V1_Adult_Mouse_Brain.h5ad',
        'image_path': 'tests/data/spatial/img.jpg',
        'annot_key': 'cluster',
    }}
fs_cfg = {'binspect_kmeans': [2000]}
cl_cfg = {'spaGCN': 1}

# , 'binspect_rank': [2000], 'spatialDE': [1000], 'SPARKX': [500]
run_bench(data_cfg, fs_cfg, cl_cfg, ['ARI', 'NMI'], modality='spatial', clean_cache=True)
# rm_cache("./cache")