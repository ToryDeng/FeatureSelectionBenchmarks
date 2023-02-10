# -*- coding: utf-8 -*-
# @Time : 2023/2/6 16:55
# @Author : Tory Deng
# @File : _load.py
# @Software: PyCharm
import os
from math import e
from typing import Dict, Literal
from typing import Union

import anndata as ad
import cv2
import scanpy as sc
from loguru import logger

from ._io import read_h5ad, write_adata_as_cache, read_adata_from_cache
from ._preprocess import make_unique, prefilter_special_genes, clean_var_names, clean_annotations
from .._utils import console


def _load_adata(data_name: str, data_props: Dict[str, Union[os.PathLike, str]], modality: Literal['scrna', 'spatial'], preprocess: bool = True,):
    logger.info("Reading adata...")
    adata = read_adata_from_cache(data_name)
    if isinstance(adata, ad.AnnData):
        console.print(f"Using cached adata and skip preprocessing. "
                      f"Shape: [yellow]{adata.n_obs}[/yellow] cells, [yellow]{adata.n_vars}[/yellow] genes.")
        return adata
    elif adata is None:
        logger.opt(colors=True).info(
            f"No cache for <magenta>{data_name}</magenta>. Trying to read h5ad from the given path in config..."
        )
        assert 'adata_path' in data_props.keys(), "Not found path to adata."
        adata = read_h5ad(data_props['adata_path'])
        if preprocess:
            console.print(f"Before QC: [yellow]{adata.n_obs}[/yellow] cells and [yellow]{adata.n_vars}[/yellow] genes.")
            clean_var_names(adata)
            make_unique(adata)
            if 'annot_key' in data_props.keys():
                to_remove = data_props['to_remove'] if 'to_remove' in data_props.keys() else None
                to_replace = data_props['to_replace'] if 'to_replace' in data_props.keys() else None
                clean_annotations(adata, data_props['annot_key'], to_remove, to_replace)
            # quality control
            prefilter_special_genes(adata)
            sc.pp.filter_genes(adata, min_cells=10)
            if modality == 'scrna':  # do not filter spots when 'modality' == 'spatial'
                sc.pp.filter_cells(adata, min_genes=200)
            console.print(f"After QC: [yellow]{adata.n_obs}[/yellow] cells and [yellow]{adata.n_vars}[/yellow] genes.")
            # store to adata.raw
            adata.raw = adata
            # normalization
            sc.pp.normalize_per_cell(adata)
            adata.layers['normalized'] = adata.X.copy()
            sc.pp.log1p(adata, base=e)
            # store data name
            adata.uns['data_name'] = data_name
            # store annot_key and batch_key
            adata.uns['annot_key'] = None if 'annot_key' not in data_props.keys() else data_props['annot_key']
            adata.uns['batch_key'] = None if 'batch_key' not in data_props.keys() else data_props['batch_key']
            # store the spot shape for spatial transcriptomics
            if 'shape' in data_props.keys():
                adata.uns['shape'] = data_props['shape']
            # save the cache
            write_adata_as_cache(adata, data_name)
        else:
            console.print(f"Skip preprocessing [yellow]{adata.n_obs}[/yellow] cells and [yellow]{adata.n_vars}[/yellow] genes.")
        return adata


def load_data(data_name: str, data_props: Dict[str, Union[os.PathLike, str]], modality: Literal['scrna', 'spatial'], preprocess: bool = True,):
    console.rule('[bold red]' + data_name)
    adata = _load_adata(data_name, data_props, modality, preprocess)
    # read image
    if modality == 'spatial':
        if 'image_path' in data_props.keys():
            img = cv2.imread(data_props['image_path'])
            logger.info("Image has been loaded.")
            return adata, img
        else:
            logger.info("'image_path' is not given. Image data not found.")
            return adata, None
    else:
        if 'image_path' in data_props.keys():
            logger.warning("The image will not be loaded when using scRNA-seq data.")
        return adata, None

