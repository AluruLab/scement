#!/usr/bin/env Rscript
rargs <- commandArgs()
args <- commandArgs(trailingOnly = TRUE)
rctypes <- c("50K", "350K", "500K", "700K", "full")
out_dir <- "fastinteg_out/"
seurat_objdir <- "/project/spc/i3/covid-scrnaseq/h5ad_full/"
if (length(args) < 2) {
  cat("Run FastIntegration for COVID data.", "\n")
  cat("Usage ", rargs, " < NCELLS > < NCORES >", "\n")
  cat("< NCELLS > should be one of : ", rctypes, "\n")
  stop("At least two arguments must be supplied (ncells, ncores)\n",
       call. = FALSE)
}
cat("Running with args : ", args, "\n")
ndata <- args[1]
ncores <- as.integer(args[2])
cat("Setting dataset run : ", ndata, ncores, "\n")
tmp_dir <- paste0(out_dir, "/", ndata)
cat("IN : ", seurat_objdir, ndata, "OUT : ", out_dir, tmp_dir, "\n")


library(Matrix)
library(ggplot2)
library(scales)
library(Seurat)
library(plotly)
library(grid)
library(tidyverse)
library(dplyr)
library(patchwork)
library(pbmcapply)
library(FastIntegration)


fast_integrate <- function(sobjects_list, tmp_dir, ftmp_dir, ncores) {
  BuildIntegrationFile(rna.list = sobjects_list, tmp.dir = tmp_dir,
                       nCores = ncores)
  FastFindAnchors(tmp.dir = tmp_dir, nCores = ncores)
  #
  genes <- readRDS(paste0(ftmp_dir, "raw/1.rds"))
  genes <- rownames(genes)
  gseq <- seq_len(length(genes))
  idx <- split(gseq, cut(gseq, ncores, labels = FALSE))
  pbmclapply(1:ncores, function(i) {
    rna_integrated <- FastIntegration(tmp.dir = tmp_dir, npcs = 1:30,
                                      slot = "data",
                                      features.to.integrate = genes[idx[[i]]])
    saveRDS(rna_integrated, paste0(ftmp_dir, "inte/inte_srt", i, ".rds"),
            compress = FALSE)
  }, mc.cores = ncores)
}

integrate_split_objects_big <- function(ftmp_dir, ncores) {
  # create Seurat obj with the variable features of integration
  #  (For very big dataset)
  features <- readRDS(paste0(ftmp_dir, "others/features.rds"))
  rna_data <- pbmclapply(1:ncores, function(i) {
    rna <- readRDS(paste0(ftmp_dir, "inte/inte_srt", i, ".rds"))
    rna <- rna[intersect(rownames(rna), features), ]
    return(rna)
  }, mc.cores = ncores)
  rna_data <- do.call(rbind, rna_data)
  rna_data <- CreateSeuratObject(rna_data)
  rna_data <- ScaleData(rna_data, features = features)
  rna_data <- RunPCA(rna_data, features = features, npcs = 50)
  rna_data <- FindNeighbors(rna_data, dims = 1:50)
  rna_data <- FindClusters(rna_data, graph.name = "RNA_snn", algorithm = 2)
  rna_data <- RunUMAP(rna_data, dims = 1:50)
  rna_data
}

integrate_split_objects <- function(ftmp_dir, ncores) {
  # Select varibale gene based on integrated data
  #   (For dataset with less than 100 samples)
  features <- readRDS(paste0(ftmp_dir, "others/features.rds"))
  rna_data <- pbmclapply(1:ncores, function(i) {
    rna <- readRDS(paste0(ftmp_dir, "inte/inte_Healthy", i, ".rds"))
    return(rna)
  }, mc.cores = 4)

  rna_data <- do.call(rbind, rna_data)
  rna_data <- CreateSeuratObject(rna_data)
  rna_data <- FindVariableFeatures(rna_data, nfeatures = 2000)
  features <- VariableFeatures(rna_data)
  rna_data <- ScaleData(rna_data, features = features)
  rna_data <- RunPCA(rna_data, features = features)
  rna_data <- FindNeighbors(rna_data, dims = 1:50)
  rna_data <- FindClusters(rna_data, resolution = 0.5, algorithm = 2)
  rna_data <- RunUMAP(rna_data, dims = 1:50)
  rna_data
}

fastinteg_wflow <- function(sobjects_list, tmp_dir, ncores) {
  ftmp_dir <- paste0(tmp_dir, "/FastIntegrationTmp/")
  #
  # sobjects_list < normalize_seurat_objects(sobjects_list)
  #
  start_time <- Sys.time()
  fast_integrate(sobjects_list, tmp_dir, ftmp_dir, ncores)
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat("FAST INTEGRATION PROC TIME", time_taken, units(time_taken),  "\n")
  #
  #
  start_time <- Sys.time()
  ncells <- Reduce(sum, lapply(sobjects_list,
                               function(x) length(Cells(x))))
  rna_data <- if (ncells > 50000) {
    integrate_split_objects_big(ftmp_dir, ncores)
  } else {
    integrate_split_objects(ftmp_dir, ncores)
  }
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat("COMBINE SPLIT TIME", time_taken, units(time_taken),  "\n")
  saveRDS(rna_data, paste0(ftmp_dir, "inte/rna.integrated.final.rds"))
}

read_seurat_object <- function(rds_file, i) {
  sdx <- readRDS(rds_file)
  new_names <- paste0(Cells(sdx), "--", i)
  sdx %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    RenameCells(new.names=new_names)
}

start_time <- Sys.time()

if (ndata == "10K") {
  use_sample <- c("S-S062", "S-M010-3", "S-S054")
}

if (ndata == "25K") {
  use_sample <-  c("S-S035-1", "S-S037", "S-S013", "S-M077", "S-S036-3")
}


if (ndata == "50K") {
  use_sample <- c("S-M018", "S-S050", "S-M066", "S-S013", "S-S035-1",
                  "S-S001-2", "S-S081", "S-M037")
}

if (ndata == "150K") {
  use_sample <- c("S-M061-2", "S-M041-1", "S-S021-4", "S-M044-1", "S-M043-2",
                 "S-M053", "S-S018", "S-S036-3", "S-S052", "S-S035-3",
                 "S-M056", "S-S087-2", "S-M049", "S-M020", "S-M001",
                 "S-S082", "S-M035-1", "S-M012", "S-S083", "S-S050",
                 "S-S027", "S-M018", "S-S086-2", "S-S061")
}


if (ndata == "350K") {
  use_sample <- c(
    "S-S041", "S-S049", "S-M060-1", "S-S076-1", "S-M011", "S-M010-4", "S-S080",
    "S-M051", "S-S020", "S-S013", "S-S022-2", "S-S039", "S-M018", "S-M007-2",
    "S-M027", "S-M004-6", "S-M033", "S-M014", "S-S018", "S-S026", "S-S086-2",
    "S-S031", "S-M042-2", "S-S073-1", "S-M008-2", "S-S083", "S-S021-4",
    "S-S043", "S-M010-3", "S-S077", "S-M004-3", "S-M017", "S-S021-2",
    "S-M005", "S-M004-2", "S-M058-1", "S-S036-1", "S-M056", "S-S091",
    "S-S070-2", "S-M007-4", "S-M010-2", "S-M076-2", "S-M043-1", "S-M028",
    "S-S030", "S-S001-2", "S-S023", "S-S035-4", "S-M041-2", "S-M007-6",
    "S-S021-3", "S-S085-2", "S-S046", "S-M008-1", "S-S033", "S-M040-2",
    "S-M077", "S-S056", "S-M009-6"
  )
}

if (ndata == "500K") {
  use_sample <- c(
    "S-M043-1", "S-M071", "S-M004-3", "S-M079", "S-M007-1", "S-M055",
    "S-S087-2", "S-M014", "S-M077", "S-M051", "S-S015", "S-S035-4",
    "S-S047", "S-S085-2", "S-S049", "S-M026-1", "S-S040", "S-S083",
    "S-M072", "S-M053", "S-M035-2", "S-M026-3", "S-M009-4", "S-M016",
    "S-S034", "S-S063", "S-S051", "S-S031", "S-M048", "S-S025", "S-M004-6",
    "S-S026", "S-M043-2", "S-S054", "S-S036-3", "S-M015", "S-S084", "S-M066",
    "S-S090-2", "S-S020", "S-S080", "S-S022-3", "S-S021-4", "S-M059-1",
    "S-M007-3", "S-M007-5", "S-S077", "S-M004-2", "S-M009-5", "S-M017",
    "S-M063-1", "S-M044-1", "S-M037", "S-S023", "S-S070-3", "S-S076-2",
    "S-S013", "S-M010-4", "S-M059-2", "S-S024", "S-M041-1", "S-M039-1",
    "S-S027", "S-S043", "S-M040-1", "S-M026-2", "S-S037", "S-M067",
    "S-S028", "S-M004-5", "S-M042-2", "S-M036-1", "S-M062-2", "S-M004-1",
    "S-M021", "S-M022", "S-M042-1", "S-M039-2", "S-M018", "S-M078"
  )
}

if (ndata == "700K") {
  use_sample <- c(
    "S-S022-2", "S-S044", "S-S015", "S-M020", "S-M056", "S-M036-1",
    "S-S037", "S-S013", "S-S088-2", "S-M048", "S-M007-2", "S-M068",
    "S-M025", "S-S021-5", "S-S080", "S-S090-2", "S-M007-1", "S-S032-3",
    "S-S079", "S-S046", "S-S022-4", "S-S036-2", "S-M011", "S-M043-2",
    "S-M009-6", "S-M077", "S-M004-4", "S-S050", "S-M037", "S-S089-2",
    "S-M008-1", "S-S078", "S-M023", "S-M029", "S-S043", "S-M061-2",
    "S-S057", "S-M060-1", "S-M061-1", "S-M009-3", "S-S087-2", "S-M035-2",
    "S-S066", "S-S045", "S-M010-6", "S-S091", "S-S035-2", "S-M041-1",
    "S-S021-2", "S-M063-1", "S-M021", "S-S022-3", "S-S086-2", "S-M047",
    "S-M074-2", "S-S040", "S-S027", "S-M034", "S-M067", "S-S063",
    "S-S065", "S-S036-3", "S-M018", "S-S029", "S-M001", "S-S023",
    "S-S068", "S-M007-6", "S-M022", "S-S022-5", "S-S075-1", "S-M009-4",
    "S-M059-1", "S-S069-3", "S-S026", "S-S051", "S-M007-4", "S-M015",
    "S-S021-4", "S-M004-1", "S-M010-5", "S-M039-2", "S-M013", "S-S016",
    "S-M071", "S-S048", "S-S083", "S-M009-5", "S-M031-2", "S-M042-2",
    "S-M033", "S-S014", "S-M045", "S-M066", "S-M049", "S-M007-3",
    "S-S025", "S-M009-1", "S-M040-2", "S-S067", "S-S042", "S-M054",
    "S-M006", "S-M024", "S-M007-5", "S-S054", "S-S039", "S-M030",
    "S-M038", "S-M079", "S-S074-1", "S-S028", "S-S031", "S-M064", "S-S084"
  )
}


if (ndata == "full") {
  use_sample <- c(
    "S-S070-2", "S-S070-3", "S-S069-3", "S-M056", "S-M044-1", "S-M043-1",
    "S-M048", "S-M044-2", "S-M043-2", "S-S054", "S-S056", "S-M042-1",
    "S-M041-1", "S-M049", "S-M046", "S-M047", "S-S055", "S-S057", "S-M045",
    "S-M041-2", "S-M042-2", "S-M055", "S-M053", "S-S067", "S-S065", "S-M051",
    "S-S064", "S-M054", "S-M052", "S-S068", "S-S066", "S-S059", "S-S060",
    "S-M050", "S-S061", "S-S062", "S-S063", "S-M061-1", "S-M061-2", "S-S073-1",
    "S-S074-1", "S-S074-2", "S-M062-2", "S-M058-1", "S-M058-2", "S-M063-1",
    "S-M063-2", "S-M059-1", "S-M059-2", "S-M060-1", "S-M060-2", "S-S075-1",
    "S-S076-1", "S-S076-2", "S-S035-1", "S-S035-2", "S-S035-3", "S-S035-4",
    "S-S036-1", "S-S036-2", "S-S036-3", "S-M064", "S-M066", "S-M067", "S-M068",
    "S-S091", "S-S090-2", "S-S092", "S-M076-2", "S-M074-2", "S-S088-2",
    "S-S085-2", "S-S089-2", "S-S087-2", "S-S086-2", "S-M077", "S-M078",
    "S-M079", "S-S024", "S-S026", "S-M013", "S-M014", "S-M015", "S-M016",
    "S-S023", "S-S025", "S-S027", "S-S034", "S-M025", "S-S033", "S-M028",
    "S-M029", "S-M027", "S-M026-1", "S-S032-3", "S-M026-2", "S-M026-3",
    "S-M004-1", "S-M004-2", "S-M004-3", "S-S022-1", "S-S022-2", "S-S022-3",
    "S-S022-4", "S-S022-5", "S-M010-1", "S-M010-2", "S-M010-3", "S-M010-4",
    "S-M010-5", "S-M010-6", "S-M009-1", "S-M009-2", "S-M011", "S-M012",
    "S-M005", "S-M006", "S-M004-4", "S-M004-5", "S-M004-6", "S-S013", "S-S014",
    "S-S015", "S-S016", "S-S017", "S-S018", "S-S019", "S-S020", "S-M007-1",
    "S-M007-2", "S-M007-3", "S-M007-4", "S-M007-5", "S-M007-6", "S-M008-1",
    "S-M008-2", "S-M009-3", "S-M009-4", "S-M009-5", "S-M009-6", "S-S021-1",
    "S-S021-2", "S-S021-3", "S-S021-4", "S-S021-5", "S-M018", "S-S029",
    "S-S030", "S-M023", "S-S031", "S-M024", "S-M017", "S-M019", "S-M020",
    "S-M021", "S-S028", "S-M022", "S-M001", "S-M030", "S-M031-1", "S-M031-2",
    "S-M032", "S-M033", "S-M034", "S-M035-1", "S-M035-2", "S-M036-1",
    "S-M036-2", "S-M037", "S-M038", "S-M039-1", "S-M039-2", "S-M040-1",
    "S-M040-2", "S-S001-2", "S-M069", "S-M070", "S-M071", "S-M072", "S-M073",
    "S-S077", "S-S078", "S-S079", "S-S080", "S-S081", "S-S082", "S-S083",
    "S-S084", "S-S037", "S-S038", "S-S039", "S-S040", "S-S041", "S-S042",
    "S-S043", "S-S044", "S-S045", "S-S046", "S-S047", "S-S048", "S-S049",
    "S-S050", "S-S051", "S-S052", "S-S053"
  )
}

list_filenames <- list.files(path = seurat_objdir,
                             pattern = ".rds$") %>%
  .[match(use_sample, gsub(".rds", "", .))]
cat("Input Files: ", length(list_filenames), "\n")
rna_list <- list()
for (i in seq_len(length(list_filenames))) {
  cat("Loading Dataset", i, list_filenames[i], "\n")
  rna_list[[i]] <- read_seurat_object(rds_file = paste0(seurat_objdir,
                                                        list_filenames[i]), i)
}
names(rna_list) <- use_sample
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("LOADING TIME", time_taken, units(time_taken),  "\n")

start_time <- Sys.time()
fastinteg_wflow(rna_list, tmp_dir, ncores)
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("FAST INTEG WORKFLOW TIME", time_taken, units(time_taken),  "\n")
