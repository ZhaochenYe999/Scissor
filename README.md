# Scissor: Single-Cell Identification of Subpopulations with bulk Sample phenOtype coRrelation #

### Introduction ###
`Scissor` is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data. The workflow of Scissor is shown in the following Figure:

<p align="center">
<img src=Figure_Method.jpg height="702" width="600">
</p>

### News ###
* May, 2021: Scissor version 2.1.0 is updated.  
    + Add utilities for cell level evaludations including correlation check and bootstrap (function: evaluate.cell)
* Feb, 2021: Scissor version 2.0.0 is launched.  
    + Optimize the inputs and outputs in Scissor main function
    + Add utilities for the reliability significance test (function: reliability.test)
* Jun, 2020: Scissor version 1.0.0 is launched.

### Installation ###
* Prerequisites:
Scissor is developed under R (*version >= 3.6.1*). The [Seurat](https://satijalab.org/seurat/) package (*version >= 3.2.0*) is used for loading data and preprocessing.

* Latest version: The latest developmental version of Scissor can be downloaded from GitHub and installed from source by
`devtools::install_github('sunduanchen/Scissor')`

### Manual ###
Please see https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html for details. In the R terminal, please use the command `?Scissor` to access the help documents.

### Examples ###
In our [Scissor Tutorial](https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html), we use several applications on the Lung Adenocarcinoma (LUAD) scRNA-seq cancer cells as examples to show how to execute Scissor in real applications.

### How to cite `Scissor` ###
Please cite the following manuscript:

> *Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data*. Nature Biotechnology (2021). https://doi.org/10.1038/s41587-021-01091-3.    
Duanchen Sun, Xiangnan Guan, Amy E. Moran, Ling-Yun Wu, David Z. Qian, Pepper Schedin, Mu-Shui Dai, Alexey V. Danilov, Joshi J. Alumkal, Andrew C. Adey, Paul T. Spellman and Zheng Xia<br />

### License ###
Scissor is licensed under the GNU General Public License v3.0.

Improvements and new features of Scissor will be updated on a regular basis. Please post on the [GitHub discussion page](https://github.com/sunduanchen/Scissor/discussions) with any questions.

### Modifications ###
This is a modified version of Scissor (modified by Zhaochen Ye, October 2026), distributed under the same GNU General Public License v3.0. Changed code is marked with `MOD` comments, with the original code kept as comments.

Main modifications include:
1. Preprocessing of single-cell and bulk data
2. Keeping the cell neighborhood matrix sparse
3. Parallelization of permutations in reliability test

#### Installation of this modified version ####
Prerequisites:
* R (*version >= 4.1.0*) and Seurat (*version >= 5.0.0*).
* A C++ compiler, since Scissor compiles C++ code during installation:
    + Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version.
    + macOS: Xcode command line tools (`xcode-select --install` in a terminal).
    + Linux: `gcc`/`g++` (e.g. `sudo apt install build-essential`). On an HPC cluster, load a GCC module; if R comes from conda, install the conda compilers instead (`conda install -c conda-forge compilers`).
* `preprocessCore` is a Bioconductor package, so it is installed separately before Scissor.

```r
install.packages(c("remotes", "BiocManager"))
BiocManager::install("preprocessCore")
remotes::install_github("ZhaochenYe999/Scissor")
```


New requirements for the input:
* `bulk_dataset` The bulk expression matrix (genes x samples), Scissor should be robust to choice of count/tpm/cpm due to the following quantile normalization.
* `phenotype` The bulk metadata. If `family="cox"`, it must be a dataframe with sample names as the rownames and the columns being survival time and censoring status, named as `time` and `status`; if `family="gaussian" or "binomial"`, it must be a named vector, matching the input sample names.
* `sc_dataset` Ideally a Seurat (v5) object. Scissor will rerun clustering for the first run, including HVG identification, batch correction, and cell-cell graph construction.