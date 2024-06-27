library(dplyr)
library(Seurat)
library(scater)
library(patchwork)
library(batchelor)
library(scran)
library(scuttle)
library(anndata)
library(zellkonverter)
library(SeuratDisk)
#library(tictoc)
source("./patch_seuratdisk.R")

# Generate analyses with fastmnn
# source data is from DATA_DIR, and output is ANALYSIS_RESULTS_DIR
DATA_DIR <- "../../data/aortic_valve"
ANALYSIS_RESULTS_DIR <- "./avalve_out/"

transpose_dgRMatrix <- function(inmat) {
    if (class(inmat) != "dgRMatrix")
        stop("inmat is not of class dgRMatrix")
    out <- new("dgCMatrix", i = inmat@j, p = inmat@p,
               x = inmat@x, Dim = rev(inmat@Dim),
               Dimnames = rev(inmat@Dimnames))
    out
}

run_fastmnn_avalve_dset_integration <- function(data_dir, out_dir) {
    #tic("PREP DATASET TIME")
    start.time <- Sys.time()
    batches <- c("H1", "H2", "C3", "C4")
    avalve_rc_sces <- lapply(batches, function(bx) {
                        h5ad_file <- paste(data_dir,
                                          paste(bx, ".h5ad", sep = ""),
                                          sep = "/")
                        ad_object <- read_h5ad(h5ad_file)
                        transpose_dgRMatrix(ad_object$X) %>%
                        CreateSeuratObject(project = bx, min.cells = 3,
                                           min.features = 200) %>%
                        as.SingleCellExperiment()
               })
    names(avalve_rc_sces) <- batches

    for (i in seq_along(avalve_rc_sces)) {
      sce_obj <- avalve_rc_sces[[i]]
      sce_feat <- rownames(sce_obj)
      #sce_obj <- sce_obj[-sort(match(pp_genes, cpt_rc_features)), ]
      #mito_features <- c(grep("ATMG", sce_feat), grep("ATCG", cpt_rc_features))
      #mito_features <- c(grep("ATMG", sce_feat), grep("ATCG", cpt_rc_features))
      mito_features <- grep("mt-", sce_feat)
      sce_obj <- addPerCellQC(sce_obj, subsets = list(Mito = mito_features))
      sce_qc <- quickPerCellQC(colData(sce_obj),
                               sub.fields = "subsets_Mito_percent")
      sce_obj <- sce_obj[, !sce_qc$discard]
      avalve_rc_sces[[i]] <- logNormCounts(sce_obj)
    }
    common_genes <- Reduce(intersect, lapply(avalve_rc_sces, rownames))

    for (i in seq_along(avalve_rc_sces)) {
      sce_obj <- avalve_rc_sces[[i]]
      avalve_rc_sces[[i]] <- sce_obj[common_genes, ]
    }
    print("Constructed Objects.")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("PREP INTEGRATION TIME", time.taken, units(time.taken),  "\n")
    start.time <- Sys.time()
    #toc()
    #tic("INTEGRATION TIME")
    mbat_out <- do.call(multiBatchNorm, avalve_rc_sces)
    model_gene_vars <- lapply(mbat_out, function(sce_obj) {
                               modelGeneVar(sce_obj)
                               })
    combined_gene_vars <- do.call(combineVar, model_gene_vars)
    chosen_hvgs <- getTopHVGs(combined_gene_vars, n = 25000)
    print("Generate HVGs.")
    combined_sce <- do.call(correctExperiments,
                             append(mbat_out, list(PARAM = NoCorrectParam())))
    #combined_sce <- runPCA(combined_sce, subset_row=chosen_hvgs)
    #combined_sce <- runUMAP(combined_sce, dimred="PCA")
    #combined_sce <- runTSNE(combined_sce, dimred="PCA")
    print("Completed correctExperiments")
    fastmnn_out <- do.call(fastMNN,
                             append(mbat_out, list(subset.row = chosen_hvgs)))
    fastmnn_out <- runUMAP(fastmnn_out, dimred = "corrected")

    print("Completed FastMNN.")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("INTEGRATION TIME", time.taken, units(time.taken),  "\n")
    start.time <- Sys.time()
    #toc()
    #tic("POST-INTEGRATION TIME")
    h5srt_file <- paste(out_dir, "fastmnn_avalve.h5Seurat", sep = "/")
    hd5out_file <- paste(out_dir, "fastmnn_avalve.h5ad", sep = "/")
    reducedDims(combined_sce) <- reducedDims(fastmnn_out)
    reducedDimNames(combined_sce) <- c("corrected", "FUMAP")
    fmnn_srt <- as.Seurat(combined_sce)
    SaveH5Seurat(fmnn_srt, h5srt_file, overwrite = TRUE)
    Convert(h5srt_file, "h5ad", overwrite = TRUE)
    #fmnn_adata <- zellkonverter::SCE2AnnData(fastmnn_out)
    #anndata::write_h5ad(fmnn_adata, filename = hd5out_file)

    rdsout_file <- paste(out_dir, "fastmnn_avalve.rds", sep = "/")
    saveRDS(fastmnn_out, rdsout_file)
    print("Saved FastMNN objects.")
    #toc()
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("POST-INTEGRATION TIME", time.taken, units(time.taken), "\n")
}

options(echo = TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(args[1] == "avalve") {
  run_fastmnn_avalve_dset_integration(DATA_DIR, ANALYSIS_RESULTS_DIR)
}
