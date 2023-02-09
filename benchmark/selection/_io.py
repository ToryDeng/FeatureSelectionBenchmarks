# -*- coding: utf-8 -*-
# @Time : 2023/2/6 19:45
# @Author : Tory Deng
# @File : _io.py
# @Software: PyCharm
import os

import numpy as np
from loguru import logger


@logger.catch
def write_genes_as_cache(genes: np.ndarray, data_name: str, fs_method: str):
    genes_dir = f"./cache/selected_genes/{data_name}/{fs_method}/"
    if not os.path.exists(genes_dir):
        os.makedirs(genes_dir)
    np.save(os.path.join(genes_dir, f"{genes.shape[0]}.npy"), genes, allow_pickle=True)
    logger.info("Selected genes have been cached.")


@logger.catch
def read_genes_from_cache(data_name: str, fs_method: str, n_genes: int):
    genes_path = f"./cache/selected_genes/{data_name}/{fs_method}/{n_genes}.npy"
    if os.path.exists(genes_path):
        logger.info("Loading cached genes...")
        return np.load(genes_path, allow_pickle=True)
    else:
        return None

