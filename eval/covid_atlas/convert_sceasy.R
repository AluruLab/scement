#!/usr/bin/env Rscript
rargs <- commandArgs()
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  cat("Rscript to convert scanpy H5AD to Seurat RDS using sceasy.", "\n")
  cat("Usage : Rscript", rargs, " < H5AD_FILE >", "\n")
  cat("Ouput is placed in < FILE BASENAME >.rds.", "\n")
  stop("At least one argument must be supplied (h5ad file)\n", call. = FALSE)
}
library(stringr)
#
input_file <- args[1]
rds_file <- str_replace(input_file, ".h5ad", ".rds")
cat("IN : ", input_file, "OUT : ", rds_file, "\n")
library(sceasy)
library(reticulate)
library(Seurat)
library(Azimuth)
#
# use_condaenv("scanpy")
loompy <- reticulate::import("loompy")
#options(future.globals.maxSize=891289600)
options(future.globals.maxSize=1073741824)
sceasy_convert(input_file, rds_file)
sceasy::convertFormat(input_file, from = "anndata", to = "seurat",
                      outFile = rds_file)
