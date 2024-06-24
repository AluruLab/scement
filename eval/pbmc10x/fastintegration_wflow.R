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

build_seurat_objects <- function(input_matrices, project_names) {
  output_objects <- vector("list", length(input_matrices))
  for (i in seq_along(input_matrices)){
    in_data <- Read10X(input_matrices[[i]])  %>%
      CreateSeuratObject(project = project_names[i], min.cells = 3,
                         min.features = 200)
    output_objects[[i]] <- in_data
  }
  output_objects
}

normalize_seurat_objects <- function(sobjects_list) {
  overlapped_genes <- Reduce(intersect, lapply(sobjects_list, rownames))
  for (i in seq_along(sobjects_list)) {
    sobjects_list[[i]] <- subset(sobjects_list[[i]],
                                 features = overlapped_genes)
    sobjects_list[[i]] <- NormalizeData(sobjects_list[[i]])
    sobjects_list[[i]] <- FindVariableFeatures(sobjects_list[[i]])
    new_names <- paste0(Cells(sobjects_list[[i]]), "--", i)
    sobjects_list[[i]] <- RenameCells(sobjects_list[[i]],
                                      new.names = new_names)
  }
  sobjects_list
}

fast_integrate <- function(sobjects_list, tmp_dir, ftmp_dir, ncores) {

  BuildIntegrationFile(rna.list = sobjects_list, tmp.dir = tmp_dir,
                       nCores = ncores)
  FastFindAnchors(tmp.dir = tmp_dir, nCores = ncores)

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
  #
  sobjects_list < normalize_seurat_objects(sobjects_list)
  #
  #
  fast_integrate(sobjects_list, tmp_dir, ftmp_dir, ncores)
  #
  #
  ncells <- Reduce(sum, lapply(sobjects_list, function(x) length(Cells(x))))
  rna_data <- if (ncells > 50000) {
    integrate_split_objects_big(ftmp_dir, ncores)
  } else {
    integrate_split_objects(ftmp_dir, ncores)
  }
  saveRDS(rna_data, paste0(ftmp_dir, "inte/rna.integrated.final.rds"))
}
