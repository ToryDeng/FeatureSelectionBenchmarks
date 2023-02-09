# -*- coding: utf-8 -*-
# @Time : 2023/2/6 0:21
# @Author : Tory Deng
# @File : _preprocess.py
# @Software: PyCharm
import re
from typing import List, Dict, Optional

import anndata as ad
import numpy as np
from loguru import logger


@logger.catch
def clean_var_names(adata: ad.AnnData):
    adata.var['original_name'] = adata.var_names
    logger.debug("The original variable names have been saved to `adata.var['original_name']`.")
    gene_names = adata.var_names.to_numpy()
    regex = re.compile(pattern='[-_:+()|]')
    vreplace = np.vectorize(lambda x: regex.sub('.', x), otypes=[str])
    adata.var_names = vreplace(gene_names)


@logger.catch
def make_unique(adata: ad.AnnData):
    if adata.obs_names.has_duplicates:
        logger.debug("Observation names have duplicates. Making them unique...")
        adata.obs_names_make_unique(join='.')
    if adata.var_names.has_duplicates:
        logger.debug("Variables names have duplicates. Making them unique...")
        adata.var_names_make_unique(join='.')
    logger.info("Observation names and Variables names are all unique now.")


@logger.catch
def clean_annotations(adata: ad.AnnData, annot_key: str, to_remove: Optional[List],
                      to_replace: Optional[Dict[str, str]]):
    # remove specified annotation types
    if to_remove:
        is_removed = adata.obs[annot_key].isin(to_remove)
        adata._inplace_subset_obs(~is_removed)
        logger.opt(colors=True).debug(f"Removed <yellow>{is_removed.sum()}</yellow> cells/spots.")
    # replace specified annotation types
    if to_replace:
        for replaced, to_replace in to_replace.items():
            if (n_to_replace := adata.obs[annot_key].isin(to_replace).sum()) > 0:
                adata.obs[annot_key].replace(to_replace, replaced, inplace=True)
                logger.opt(colors=True).debug(f"Replaced <yellow>{n_to_replace}</yellow> annotations with '{replaced}'.")



@logger.catch
def prefilter_special_genes(adata: ad.AnnData, Gene1Pattern: str = "ERCC", Gene2Pattern: str = "MT-"):
    drop_pattern1 = adata.var_names.str.startswith(Gene1Pattern)
    drop_pattern2 = adata.var_names.str.startswith(Gene2Pattern)
    drop_pattern = np.logical_and(~drop_pattern1, ~drop_pattern2)
    logger.opt(colors=True).info(f"Dropping <yellow>{adata.n_vars - drop_pattern.sum()}</yellow> special genes "
                                 f"from <yellow>{adata.n_vars}</yellow> genes...")
    adata._inplace_subset_var(drop_pattern)
