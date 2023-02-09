# -*- coding: utf-8 -*-
# @Time : 2023/2/7 18:49
# @Author : Tory Deng
# @File : _general.py
# @Software: PyCharm
from typing import Union, Callable, Literal, Optional

import anndata as ad
import numpy as np

from ._io import read_genes_from_cache, write_genes_as_cache
from .scrna.functions import select_features as scrna_select
from .spatial.functions import select_features as spatial_select
from .._utils import console


def generally_select_features(
        adata: ad.AnnData,
        img: Optional[np.ndarray],
        fs_method: Union[str, Callable],
        n_selected_genes: Union[int, Literal['auto']],
        modality: Literal['scrna', 'spatial'],
        random_state: int,
        **kwargs
) -> np.ndarray:
    if isinstance(fs_method, str):
        console.rule(f'FS method: [bold red]{fs_method}[/bold red], Gene Number: [bold red]{n_selected_genes}[/bold red]')
        # load cached genes or run feature selection
        selected_genes = read_genes_from_cache(adata.uns['data_name'], fs_method, n_selected_genes)
        if selected_genes is None:
            if modality == 'scrna':
                selected_genes = scrna_select(adata, fs_method, n_selected_genes, random_state)
            else:
                selected_genes = spatial_select(adata, img, fs_method, n_selected_genes, random_state)
            write_genes_as_cache(selected_genes, adata.uns['data_name'], fs_method)
    elif callable(fs_method):  # fs_method is a function
        console.rule(f'FS method: [bold red]{fs_method.__name__}[/bold red], Gene Number: [bold red]{n_selected_genes}[/bold red]')
        # load cached genes or run feature selection
        selected_genes = read_genes_from_cache(adata.uns['data_name'], fs_method.__name__, n_selected_genes)
        if selected_genes is None:
            selected_genes = fs_method(adata, n_selected_genes, **kwargs)
            write_genes_as_cache(selected_genes, adata.uns['data_name'], fs_method.__name__)
    else:
        raise NotImplementedError(f"`fs_method` should be an valid string or a function, got {type(fs_method)}.")
    return selected_genes
