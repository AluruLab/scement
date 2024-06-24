#!/usr/bin/env Rscript
rargs <- commandArgs()
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)\n", call. = FALSE)
}
library(stringr)
cat("Running with args : ", args, "\n")
input_file <- args[1]
azrds_file <- str_replace(input_file, ".h5ad", "_azimuth.rds")
mdata_file <- str_replace(input_file, ".h5ad", "_mdata.csv")
rds_file <- str_replace(input_file, ".h5ad", ".rds")
cat("Setting dataset run : ", input_file, "\n")
cat("IN : ", input_file, "OUT : ", azrds_file, mdata_file, rds_file, "\n")

library(sceasy)
library(reticulate)
library(Seurat)
library(Azimuth)
# use_condaenv("scanpy")
loompy <- reticulate::import("loompy")

sceasy_convert <- function(h5adfile, rdsfile) {
  sceasy::convertFormat(h5adfile, from = "anndata", to = "seurat",
                        outFile = rdsfile)
}

run_azimuth <- function(rds_file, azrds_file, mdata_file) {
  sobj <- readRDS(rds_file)
  sobj <- RunAzimuth(sobj, reference = "pbmcref")
  saveRDS(sobj, azrds_file)
  sobj <-  sobj@meta.data
  write.csv(sobj, mdata_file)
}

#options(future.globals.maxSize=891289600)
options(future.globals.maxSize=1073741824)
sceasy_convert(input_file, rds_file)
run_azimuth(rds_file, azrds_file, mdata_file)
