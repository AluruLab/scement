library(dplyr)
library(Seurat)
library(scater)
library(patchwork)
library(batchelor)
library(scran)
library(scuttle)
library(anndata)

args <- commandArgs(trailingOnly = TRUE)
cat("Running with args : ", args, "\n")
ndata <- as.integer(args[1])
cat("Setting No. of datasets : ", ndata, "\n")
out_dir <- "fastmnn_out"
seurat_objdir <- "../../data/FastIntegrate/Seurat_Objects/"

if (ndata == 2) {
    use_sample <- c("n06", "n07")
}

if (ndata == 3) {
    use_sample <- c("n02", "n06", "n07")
}

if (ndata == 5) {
    use_sample <- c("n01", "n02", "n06", "n07", "n08")
}

if (ndata == 9) {
    use_sample <- c("n01", "n02", "n06", "n07", "n08", "n03", "n05", "n13", "n14")
}

if (ndata == 13) {
    use_sample <- c("n01", "n02", "n06", "n07", "n08", "n03", "n05", "n09", "n10",
                    "n11", "n12", "n13", "n14")
}

if (ndata == 16) {
    use_sample <- c("n01", "n02", "n06", "n07", "n08", "n03", "n05", "n09", "n10",
                    "n11", "n12", "n13", "n14", "n70", "n71", "n74")
}

if (ndata == 17) {
    use_sample <- c("n01", "n02", "n06", "n07", "n08", "n03", "n05", "n09", "n10",
                    "n11", "n12", "n13", "n14", "n70", "n71", "n74", "n15")
}


read_seu <- function(rds_file) {
    as.SingleCellExperiment(readRDS(rds_file))
}


list_filenames <- list.files(path = seurat_objdir,
                             pattern = ".rds$") %>%
    .[match(use_sample, gsub("_seurat.rds", "", .))]
cat("Input Files: ", length(list_filenames), "\n")
rna_list <- list()
for (i in seq_len(length(list_filenames))) {
    cat("Loading Dataset", i, list_filenames[i], "\n")
    rna_list[[i]] <- read_seu(rds_file = paste0(seurat_objdir, list_filenames[i]))
}
names(rna_list) <- use_sample

start_time <- Sys.time()
for (i in seq_along(rna_list)) {
    sce_obj <- rna_list[[i]]
    sce_feat <- rownames(sce_obj)
    #sce_obj <- sce_obj[-sort(match(pp_genes, cpt_rc_features)), ]
    #mito_features <- c(grep("ATMG", sce_feat), grep("ATCG", cpt_rc_features))
    #mito_features <- c(grep("ATMG", sce_feat), grep("ATCG", cpt_rc_features))
    mito_features <- grep("mt-", sce_feat)
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
combined_sce <- do.call(correctExperiments,
                        append(mbat_out, list(PARAM = NoCorrectParam())))
# combined_sce <- runPCA(combined_sce, subset_row=chosen_hvgs)
# combined_sce <- runUMAP(combined_sce, dimred="PCA")
# combined_sce <- runTSNE(combined_sce, dimred="PCA")
# cat("Completed correctExperiments", "\n")
cat("Completed Combo : ", str(assay(combined_sce)), "\n")
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
