# -*- coding: utf-8 -*-
# @Time : 2023/2/6 22:56
# @Author : Tory Deng
# @File : _metrics.py
# @Software: PyCharm
from typing import Literal, Union

import numpy as np
import pandas as pd
from loguru import logger
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score


def compute_clustering_metrics(true_labels: Union[np.ndarray, pd.Series], pred_labels: np.ndarray, metric: Literal['ARI', 'NMI']):
    valid_label_mask = ~pd.isna(true_labels).to_numpy()  # np.isnan doesn't support strings
    if (n_nan_labels := true_labels.shape[0] - valid_label_mask.sum()) > 0:
        logger.opt(colors=True).info(
            f"Ignoring <yellow>{n_nan_labels}</yellow> NaNs in annotations during the {metric} computation."
        )
    if metric == 'ARI':
        return adjusted_rand_score(true_labels[valid_label_mask], pred_labels[valid_label_mask])
    elif metric == 'NMI':
        return normalized_mutual_info_score(true_labels[valid_label_mask], pred_labels[valid_label_mask])
    else:
        raise NotImplementedError(f"Metric {metric} is not implemented.")
