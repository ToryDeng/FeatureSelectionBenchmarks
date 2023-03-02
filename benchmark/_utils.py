# -*- coding: utf-8 -*-
# @Time : 2023/2/5 23:43
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import os
import shutil
import sys
from typing import Union, Literal, Optional

from loguru import logger
from rich.console import Console

console = Console()


class HiddenPrints:
    """
    Hide prints from terminal
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr


def rm_cache(path: Union[os.PathLike, str]):
    if os.path.exists(path):
        shutil.rmtree(path)
        logger.info(f"{path} has been deleted.")
    else:
        logger.warning(f"{path} not found. Skip deleting the directory.")


def set_logger(verbosity: Literal[0, 1, 2] = 1, log_path: Optional[Union[os.PathLike, str]] = None):
    """
    Set the verbosity level.

    Parameters
    ----------
    verbosity
        0: only print warnings and errors
        1: also print info
        2: also print debug messages
    log_path
        Path to the log file
    """
    def formatter(record: dict):
        if record['level'].name in ('DEBUG', 'INFO'):
            return "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | " \
                   "<level>{level: <5}</level> | " \
                   "<level>{message}\n</level>"
        else:
            return "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | " \
                   "<level>{level: <8}</level> | " \
                   "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}\n</level>"

    level_dict = {0: 'WARNING', 1: 'INFO', 2: 'DEBUG'}
    logger.remove()
    if log_path:
        logger.add(log_path, encoding='utf-8')
    else:
        logger.add(sys.stdout, colorize=True, level=level_dict[verbosity], format=formatter)

