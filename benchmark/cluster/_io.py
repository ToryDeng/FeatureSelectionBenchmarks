# -*- coding: utf-8 -*-
# @Time : 2023/2/6 22:42
# @Author : Tory Deng
# @File : _io.py
# @Software: PyCharm
import os

import numpy as np
from loguru import logger


def write_clusters_as_cache(clusters: np.ndarray, data_name: str, fs_method: str, n_genes: int, cl_method: str, run: int):
    clusters_dir = f"./cache/clustering_result/{data_name}/{fs_method}/{n_genes}/{cl_method}/"
    if not os.path.exists(clusters_dir):
        os.makedirs(clusters_dir)
    np.save(os.path.join(clusters_dir, f"{run}.npy"), clusters, allow_pickle=True)
    logger.opt(colors=True).info(f"<magenta>{cl_method}</magenta> clustering results have been cached.")


def read_clusters_from_cache(data_name: str, fs_method: str, n_genes: int, cl_method: str, run: int):
    clusters_dir = f"./cache/clustering_result/{data_name}/{fs_method}/{n_genes}/{cl_method}/{run}.npy"
    if os.path.exists(clusters_dir):
        logger.opt(colors=True).info(f"Loading cached <magenta>{cl_method}</magenta> clustering results...")
        return np.load(clusters_dir, allow_pickle=True)
    else:
        return None
