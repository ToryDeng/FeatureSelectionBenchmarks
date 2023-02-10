# -*- coding: utf-8 -*-
# @Time : 2023/2/7 1:03
# @Author : Tory Deng
# @File : _recorder.py
# @Software: PyCharm
import os
from datetime import datetime
from itertools import product
from typing import Dict, Union, List, Literal, Callable

import numpy as np
import pandas as pd
from loguru import logger


def create_records(
        data_cfg: Dict[str, Dict[str, Union[os.PathLike, str, List, Dict[str, Union[List, str]]]]],
        fs_cfg: Dict[str, Union[List[int], List[Literal['auto']]]],
        cl_cfg: Dict[str, int],
        metrics: List[Literal['ARI', 'NMI']]
):
    col_tuples = [
        tup
        for fs_method, n_genes in fs_cfg.items()
        for tup in product((fs_method.__name__ if callable(fs_method) else fs_method,), n_genes)
    ]
    col_index = pd.MultiIndex.from_tuples(col_tuples, names=['fs_method', 'n_genes'])

    row_tuples = [
        tup
        for cl_method, n_runs in cl_cfg.items()
        for tup in product(data_cfg.keys(), (cl_method.__name__ if callable(cl_method) else cl_method,), range(n_runs))
    ]
    row_index = pd.MultiIndex.from_tuples(row_tuples, names=['dataset', 'clustering_method', 'run'])
    single_record = pd.DataFrame(
        np.full((len(row_tuples), len(col_tuples)), fill_value=np.nan, dtype=float), index=row_index, columns=col_index
    )
    return {metric: single_record.copy() for metric in metrics}


def store_metrics_to_records(
        records: Dict[str, pd.DataFrame],
        metric: Literal['ARI', 'NMI'],
        value: float,
        data_name: str,
        cl_method: Union[str, Callable],
        run: int,
        fs_method: Union[str, Callable],
        n_selected_genes: Union[int, Literal['auto']]
):
    if callable(fs_method):
        fs_method = fs_method.__name__
    if callable(cl_method):
        cl_method = cl_method.__name__
    records[metric].loc[(data_name, cl_method, run), (fs_method, n_selected_genes)] = value


def write_records(records: Dict[str, pd.DataFrame], modality: Literal['scrna', 'spatial']):
    record_name = f"{datetime.now().strftime('%Y-%m %H_%M_%S')} {modality}"
    writer = pd.ExcelWriter(f'{record_name}.xlsx')
    for metric, record in records.items():
        record.to_excel(writer, sheet_name=metric, index=True)
    writer.close()
    logger.info(f"records have been saved into './{record_name}.xlsx'.")
