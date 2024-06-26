#!/usr/bin/env Rscript
rargs <- commandArgs()
args <- commandArgs(trailingOnly = TRUE)
rctypes <- c("50K", "350K", "500K", "700K", "full")
seurat_objdir <- "../../data/covid_atlas/"
out_dir <- "fastmnn_out/"
if (length(args) == 0) {
  cat("Run FastMNN for COVID data.", "\n")
  cat("Usage ", rargs, "< NCELLS >", "\n")
  cat("< NCELLS > should be one of : ", rctypes, "\n")
  stop("At least one argument must be supplied (ncells)\n", call. = FALSE)
}
cat("Running with args : ", args, "\n")
ndata <- args[1]
cat("Setting dataset run : ", ndata, "\n")
cat("IN : ", seurat_objdir, ndata, "OUT : ", out_dir, "\n")

library(dplyr)
library(Seurat)
library(scater)
library(patchwork)
library(batchelor)
library(scran)
library(scuttle)
library(anndata)

if (ndata == "10k") {
  use_sample <- c("S-S062", "S-M010-3", "S-S054")
}

if (ndata == "25k") {
  use_sample <-  c("S-S035-1", "S-S037", "S-S013", "S-M077", "S-S036-3")
}

if (ndata == "50k") {
  use_sample <- c("S-M018", "S-S050", "S-M066", "S-S013", "S-S035-1",
                  "S-S001-2", "S-S081", "S-M037")
}

if (ndata == "150k") {
  use_sample <- c("S-M061-2", "S-M041-1", "S-S021-4", "S-M044-1", "S-M043-2",
                 "S-M053", "S-S018", "S-S036-3", "S-S052", "S-S035-3",
                 "S-M056", "S-S087-2", "S-M049", "S-M020", "S-M001",
                 "S-S082", "S-M035-1", "S-M012", "S-S083", "S-S050",
                 "S-S027", "S-M018", "S-S086-2", "S-S061")
}

if (ndata == "350k") {
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

if (ndata == "500k") {
  use_sample <- c(
    "S-M043-1", "S-M071", "S-M004-3", "S-M079", "S-M007-1", "S-M055",
    "S-S087-2", "S-M014", "S-M077", "S-M051", "S-S015", "S-S035-4", "S-S047",
    "S-S085-2", "S-S049", "S-M026-1", "S-S040", "S-S083", "S-M072", "S-M053",
    "S-M035-2", "S-M026-3", "S-M009-4", "S-M016", "S-S034", "S-S063", "S-S051",
    "S-S031", "S-M048", "S-S025", "S-M004-6", "S-S026", "S-M043-2", "S-S054",
    "S-S036-3", "S-M015", "S-S084", "S-M066", "S-S090-2", "S-S020", "S-S080",
    "S-S022-3", "S-S021-4", "S-M059-1", "S-M007-3", "S-M007-5", "S-S077",
    "S-M004-2", "S-M009-5", "S-M017", "S-M063-1", "S-M044-1", "S-M037",
    "S-S023", "S-S070-3", "S-S076-2", "S-S013", "S-M010-4", "S-M059-2",
    "S-S024", "S-M041-1", "S-M039-1", "S-S027", "S-S043", "S-M040-1",
    "S-M026-2", "S-S037", "S-M067", "S-S028", "S-M004-5", "S-M042-2",
    "S-M036-1", "S-M062-2", "S-M004-1", "S-M021", "S-M022", "S-M042-1",
    "S-M039-2", "S-M018", "S-M078"
  )
}

if (ndata == "700k") {
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

read_seu <- function(rds_file) {
  as.SingleCellExperiment(readRDS(rds_file))
}

start_time <- Sys.time()
list_filenames <- list.files(path = seurat_objdir,
                             pattern = ".rds$") %>%
  .[match(use_sample, gsub(".rds", "", .))]
cat("Input Files: ", length(list_filenames), "\n")
rna_list <- list()
for (i in seq_len(length(list_filenames))) {
  cat("Loading Dataset", i, list_filenames[i], "\n")
  rna_list[[i]] <- read_seu(rds_file = paste0(seurat_objdir, list_filenames[i]))
}
names(rna_list) <- use_sample
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("LOADING TIME", time_taken, units(time_taken),  "\n")

start_time <- Sys.time()
for (i in seq_along(rna_list)) {
  sce_obj <- rna_list[[i]]
  sce_feat <- rownames(sce_obj)
  #sce_obj <- sce_obj[-sort(match(pp_genes, cpt_rc_features)), ]
  #mito_features <- c(grep("ATMG", sce_feat), grep("ATCG", cpt_rc_features))
  #mito_features <- c(grep("ATMG", sce_feat), grep("ATCG", cpt_rc_features))
  mito_features <- grep("MT-", sce_feat)
  sce_obj <- addPerCellQC(sce_obj, subsets = list(Mito = mito_features))
  sce_qc <- quickPerCellQC(colData(sce_obj),
                           sub.fields = "subsets_Mito_percent")
  sce_obj <- sce_obj[, !sce_qc$discard]
  rna_list[[i]] <- logNormCounts(sce_obj)
}
common_genes <- Reduce(intersect, lapply(rna_list, rownames))

for (i in seq_along(rna_list)) {
  sce_obj <- rna_list[[i]]
  rna_list[[i]] <- sce_obj[common_genes, ]
}
print("Constructed Objects.")
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("PREP INTEGRATION TIME", time_taken, units(time_taken),  "\n")

start_time <- Sys.time()
mbat_out <- do.call(multiBatchNorm, rna_list)
cat("Completed MultiBatch : ", length(mbat_out), "\n")
model_gene_vars <- lapply(mbat_out, function(sce_obj) {
  modelGeneVar(sce_obj)
})
combined_gene_vars <- do.call(combineVar, model_gene_vars)
chosen_hvgs <- getTopHVGs(combined_gene_vars, n = 25000)
cat("Generate HVGs : ", length(chosen_hvgs), "\n")
# combined_sce <- do.call(correctExperiments,
#                        append(mbat_out, list(PARAM = NoCorrectParam())))
# combined_sce <- runPCA(combined_sce, subset_row=chosen_hvgs)
# combined_sce <- runUMAP(combined_sce, dimred="PCA")
# combined_sce <- runTSNE(combined_sce, dimred="PCA")
# cat("Completed correctExperiments", "\n")
# cat("Completed Combo : ", str(assay(combined_sce)), "\n")
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("COMBO TIME", time_taken, units(time_taken),  "\n")

start_time <- Sys.time()
fastmnn_out <- do.call(fastMNN, mbat_out)
# fastmnn_out <- do.call(fastMNN,
#                       append(mbat_out, list(subset.row = chosen_hvgs)))
fastmnn_out <- runUMAP(fastmnn_out, dimred = "corrected")
cat("Completed FastMNN : ", dim(fastmnn_out), str(assay(fastmnn_out)), "\n")
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("INTEGRATION TIME", time_taken, units(time_taken),  "\n")

start_time <- Sys.time()
#  start_time <- Sys.time()
#  h5srt_file <- paste(out_dir, "fastmnn_avalve.h5Seurat", sep = "/")
#  hd5out_file <- paste(out_dir, "fastmnn_avalve.h5ad", sep = "/")
#  reducedDims(combined_sce) <- reducedDims(fastmnn_out)
#  reducedDimNames(combined_sce) <- c("corrected", "FUMAP")
#  fmnn_srt <- as.Seurat(combined_sce)
#  SaveH5Seurat(fmnn_srt, h5srt_file, overwrite = TRUE)
#  Convert(h5srt_file, "h5ad", overwrite = TRUE)
#  #fmnn_adata <- zellkonverter::SCE2AnnData(fastmnn_out)
#  #anndata::write_h5ad(fmnn_adata, filename = hd5out_file)

rdsout_file <- paste0(out_dir, "/fastmnn_avalve_", ndata, ".rds")
saveRDS(fastmnn_out, rdsout_file)
cat("Saved FastMNN objects", "\n")
end_time <- Sys.time()
time_taken <- end_time - start_time
cat("POST-INTEGRATION TIME", time_taken, units(time_taken), "\n")
