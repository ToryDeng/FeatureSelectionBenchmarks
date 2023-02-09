# -*- coding: utf-8 -*-
# @Time : 2023/2/5 23:43
# @Author : Tory Deng
# @File : _utils.py
# @Software: PyCharm
import os
import shutil
import sys
from typing import Union

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

