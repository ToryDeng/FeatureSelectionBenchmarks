# FeatureSelectionBenchMarks

This is the code repository of the feature (gene) selection benchmark in both scRNA-seq and spatial transcriptomics.

## Software features

After a simple configuration, you can run the benchmark (including data loading, quality control, feature selection, and cell clustering/domain detection) in **one single line of code**:

```python
from benchmark.run_benchmark import run_bench


# configure the dataset information
data_cfg = {
    'your_data_name': {
        'adata_path': 'path/to/h5ad/file',
        'annot_key': 'annotation_name',
    }}
# configure feature selection methods and numbers of selected features
fs_cfg = {'feature_selection_method': [1000, 2000]}
# configure clustering methods and numbers of runs
cl_cfg = {'clustering_method': 2}
# run the benchmark in one line of code
run_bench(data_cfg, fs_cfg, cl_cfg, modality='scrna', metrics=['ARI', 'NMI'])
```

The evaluation results will be automatically saved as an XLSX file in the working directory with name like this:

```text
2023-02 14_54_32 scrna.xlsx
```

Other software features are:

- Automatically save the results of each step (preprocessed data, selected features, and cluster labels)
- Reload the cached genes and cluster labels when you use the same data (specified by the data name)
- Support custom feature selection and cell clustering/domain detection methods
- Present detailed and pretty logging messages based on [rich](https://github.com/Textualize/rich) and [loguru](https://github.com/Delgan/loguru) (see examples in tutorial)

## Currently supported methods

### scRNA-seq

#### Feature selection

| Name  | Language | Reference |
| :---: | :---:    | :---:     |
| GeneClust | Python | [paper](https://doi.org/10.1093/bib/bbad042)
| vst   | Python   | [paper](https://doi.org/10.1016/j.cell.2019.05.031) |
| mvp   | Python   | [paper](https://www.nature.com/articles/nbt.3192) |
| triku | Python   | [paper](https://doi.org/10.1093/gigascience/giac017) |
| GiniClust3|Python| [paper](https://doi.org/10.1186/s12859-020-3482-1) |
| SC3   | Python        | [paper](https://doi.org/10.1038/nmeth.4236) |
| scran | R        | [paper](https://doi.org/10.1186/s13059-016-0947-7) |
| FEAST | R        | [paper](https://doi.org/10.1093/bib/bbab034) |
| M3Drop | R       | [paper](https://doi.org/10.1093/bioinformatics/bty1044) |
| scmap  | R        | [paper](https://doi.org/10.1038/nmeth.4644) |
| deviance    | R        | [paper](https://doi.org/10.1186/s13059-019-1861-6) |
| FEAST       | R        | [paper](https://doi.org/10.1093/bib/bbab034) |
| sctransform | R        | [paper](https://doi.org/10.1186/s13059-019-1874-1) |

#### Cell clustering

| Name  | Language | Reference |
| :---: | :---:    | :---:     |
| SC3s | Python  | [paper](https://doi.org/10.1186/s12859-022-05085-z) |
| Seurat | R       | [paper](https://doi.org/10.1016/j.cell.2021.04.048) |
| SHARP  | R       | [paper](http://www.genome.org/cgi/doi/10.1101/gr.254557.119) |
| TSCAN | R       | [paper](https://doi.org/10.1093/nar/gkw430) |
| CIDR | R       | [paper](https://doi.org/10.1186/s13059-017-1188-0) |


### Spatial transcriptomics
#### Feature selection

| Name  | Language | Reference |
| :---: | :---:    | :---:     |
| SpatialDE |Python| [paper](https://doi.org/10.1038/nmeth.4636v) |
| SPARK-X   |   R  | [paper](https://doi.org/10.1186/s13059-021-02404-0) |
| Giotto    |   R  | [paper](https://doi.org/10.1186/s13059-021-02286-2)

#### Domain detection

| Name  | Language | Reference |
| :---: | :---:    | :---:     |
| SpaGCN | Python  | [paper](https://doi.org/10.1038/s41592-021-01255-8) |
| stLearn | Python  | [paper](https://doi.org/10.1101/2020.05.31.125658) |
| STAGATE | Python  | [paper](https://doi.org/10.1038/s41467-022-29439-6) |

## Requirements

### R packages

This benchmark is written in Python and calls R functions through `rpy2`. If you want to use some methods implemented with R language, please install the corresponding R packages.

### Python packages

- anndata>=0.8.0
- numpy>=1.21.6
- setuptools>=59.5.0
- anndata2ri>=1.1
- sc3s>=0.1.1
- scanpy>=1.9.1
- loguru>=0.6.0
- rpy2>=3.5.6
- sklearn>=0.0.post2
- scikit-learn>=1.2.0
- SpaGCN>=1.2.5
- torch>=1.13.1
- stlearn>=0.4.11
- pandas>=1.5.2
- opencv-python>=4.6.0
- scipy>=1.9.3
- rich>=13.0.0
- triku>=2.1.4
- statsmodels>=0.13.5
- SpatialDE>=1.1.3
- STAGATE_pyG>=1.0.0

## Installation
```commandline
git clone https://github.com/ToryDeng/FeatureSelectionBenchmarks
cd FeatureSelectionBenchmarks/
python3 setup.py install --user
```

## Tutorial
- The tutorial about how to run the benchmarks: [tutorials/run_benchmarks.ipynb](https://github.com/ToryDeng/FeatureSelectionBenchmarks/blob/main/tutorials/run_benchmarks.ipynb)
- The tutorial about how to read the records: [tutorials/read_records.ipynb](https://github.com/ToryDeng/FeatureSelectionBenchmarks/blob/main/tutorials/read_records.ipynb)



