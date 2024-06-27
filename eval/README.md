# Evaluation of SCEMENT and other integration methods

This directory contains the scripts to run the methods run for
evaluation of SCEMENT provided in the paper. 
In the paper, we used human aortic valve and A. thaliana to evaluate 
the quality of the methods, while we used pbmc datasets from COVID-19 
atlas and public databases (including 10X datasest) for evaluation of 
scalability. 

The scripts for each of these evaluation is given in folders avalve,
athaliana, covid_atlas and pbmc10x respectively. Since the software
used for evaluation are written in either R or python, the 
 evaluation scripts are written either in python or R.

## Data

The raw data is pre-processed and converted into anndata h5ad files.
The pre-processed data can be downloaded from the zenodo resource [https://zenodo.org/doi/10.5281/zenodo.11521687](https://zenodo.org/doi/10.5281/zenodo.11521687)  as decribed in data/README file.
Evaluation scripts assume that the data is downloaded and unzipped in the data directory before running R/python scripts. 

## Dependencies

The evaluation scripts require that the following R and python libararies need to be installed.

### R dependencies

Evaluation scripts depend on the following R packages. These packages can be installed from the Bioconductor/CRAN repositories.

   - Azimuth
   - BPCells
   - FastIntegration
   - Matrix
   - Seurat
   - SeuratDisk
   - SeuratObject
   - anndata
   - batchelor
   - dplyr
   - ggplot2
   - grid
   - patchwork
   - pbmcapply
   - plotly
   - reticulate
   - scales
   - scater
   - sceasy
   - scran
   - scuttle
   - stringr
   - tictoc
   - tidyverse
   - zellkonverter

### python dependencies

Evaluation scripts depend on the following python libraries, which can be installed using pip. Creating a virtual environment using conda/virtualenv is recommended.

   - anndata
   - numpy
   - pandas
   - scanorama
   - scanpy
   - scement (Library at this source)
   - scement
   - scib
   - scipy

## Outputs

Outputs of the evaluation are generated within the respective directories - avalve, athaliana, covid_atlas and pbmc10x.
