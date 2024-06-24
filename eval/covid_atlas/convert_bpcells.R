#!/usr/bin/env Rscript
rargs <- commandArgs()
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  cat("Rscript to convert Seurat RDS file to BPB Cells data.", "\n")
  cat("Usage : ", rargs, "< RDS_FILE > ", "\n")
  cat("Ouput is placed in < FILE BASENAME >_data_BP directory", "\n")
  stop("At least one argument must be supplied (rds file)\n", call. = FALSE)
}
#
library(stringr)
input_file <- args[1]
out_bp_dir <- str_replace(input_file, ".rds", "_data_BP")
cat("IN : ", input_file, "OUT : ", out_bp_dir, "\n")
#
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)

# options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1073741824)
rdx <- readRDS(input_file)
dgmat <- rdx@assays[["RNA"]]@counts
write_matrix_dir(mat = dgmat, dir = out_bp_dir, overwrite = TRUE)
