# -*- coding: utf-8 -*-
# @Time : 2023/2/9 20:50
# @Author : Tory Deng
# @File : setup.py
# @Software: PyCharm
from setuptools import find_packages
from setuptools import setup

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setup(
    name='FeatureSelectionBenchmarks',
    version='1.0.0',
    description='Benchmarking feature selection methods for scRNA-seq and spatially resolved transcriptomics',
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    author='Tao Deng',
    author_email='taodeng@link.cuhk.edu.cn',
    url='https://github.com/ToryDeng/FeatureSelectionBenchmarks',
    license='GPL v3',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    python_requires=">=3.9",
    install_requires=[
        'anndata>=0.8.0',
        'numpy>=1.21.6',
        'scanpy>=1.9.1',
        'loguru>=0.6.0',
        'anndata2ri>=1.1',
        'sc3s>=0.1.1',
        'rpy2>=3.5.6',
        'scikit-learn>=1.2.0',
        'SpaGCN>=1.2.5',
        'torch>=1.13.1',
        'stlearn>=0.4.11',
        'pandas>=1.5.2',
        'opencv-python>=4.6.0',
        'scipy>=1.9.3',
        'rich>=13.0.0',
        'triku>=2.1.4',
        'statsmodels>=0.13.5',
        'SpatialDE>=1.1.3'
    ],
    packages=find_packages(exclude=('tests', 'figures', 'data', 'docs', 'notebooks', 'tutorials')),
    zip_safe=False,
)
