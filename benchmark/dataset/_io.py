# -*- coding: utf-8 -*-
# @Time : 2023/2/5 23:00
# @Author : Tory Deng
# @File : _io.py
# @Software: PyCharm
import os
from typing import Union

import anndata as ad
import scanpy as sc
from loguru import logger

from ._utils import is_normalized, to_dense


@logger.catch
def read_h5ad(path: Union[os.PathLike, str], use_raw: bool = True) -> ad.AnnData:
    adata = sc.read_h5ad(path)
    if use_raw and adata.raw is not None:
        logger.info("Using raw data in `adata.raw`...")
        adata = adata.raw.to_adata()
        to_dense(adata)
        if is_normalized(adata):
            raise ValueError("Found normalized data in `adata.raw`.")
    else:
        to_dense(adata)
    return adata


@logger.catch
def write_adata_as_cache(adata:ad.AnnData, data_name: str):
    if not os.path.exists("./cache/h5ad/"):
        os.makedirs("./cache/h5ad/")
    adata.write_h5ad(f"./cache/h5ad/{data_name}.h5ad")
    logger.debug(f"adata has been written to ./cache/h5ad/{data_name}.h5ad")


@logger.catch
def read_adata_from_cache(data_name: str):
    if os.path.exists(f"./cache/h5ad/{data_name}.h5ad"):
        adata = read_h5ad(f"./cache/h5ad/{data_name}.h5ad", use_raw=False)
        # if (n_str_nan := (adata.obs['spatialLIBD'] == 'nan').sum()) > 0:
        #     logger.opt(colors=True).info(
        #         f"Found <yellow>{n_str_nan}</yellow> NaNs in annotations with string format ('nan'). "
        #         f"Converting them to np.nan...")
        #     adata.obs['spatialLIBD'].replace('nan', np.nan, inplace=True)
        logger.opt(colors=True).info(f"Loaded cahed <magenta>{data_name}</magenta> dataset.")

        return adata
    else:
        return None

