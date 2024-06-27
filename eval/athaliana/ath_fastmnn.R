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
source("./patch_seuratdisk.R")

DATA_DIR <- "../../data/athaliana/"
PROTPLAST_GENE_LIST <-  "../../data/meta/Protoplasting_DEgene_FC2_list.txt"
ANAYSIS_OUT_DIR <- "./athaliana_out/"

transpose_dgRMatrix <- function(inmat) {
    if (class(inmat) != "dgRMatrix")
        stop("inmat is not of class dgRMatrix")
    out <- new("dgCMatrix", i = inmat@j, p = inmat@p,
               x = inmat@x, Dim = rev(inmat@Dim),
               Dimnames = rev(inmat@Dimnames))
    out
}


run_gala_dset <- function(data_dir, proto_file, out_dir) {
    gala_file <- paste(data_dir, "GSE158761/GSE158761.h5ad", sep="/")
    gala_ds <- read_h5ad(gala_file)
    gala_ds_matrix <- gala_ds$X
    gala_ds_obs <- gala_ds$obs
    gala_ds_counts <- transpose_dgRMatrix(gala_ds_matrix)
    gala_ds_batches <- levels(gala_ds_obs[, "batch"])
    gala_list <- lapply(gala_ds_batches, function(bx) {
        gala_ds_bcodes <- gala_ds_obs[gala_ds_obs["batch"] == bx, "Assay"]
        gala_ds_bcounts <-  gala_ds_counts[, gala_ds_bcodes]
        gala_ds_bcounts %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            as.SingleCellExperiment()
    })
    names(gala_list) <- gala_ds_batches
    saveRDS(gala_list, file=paste(out_dir, "fmnn_gala_list.rds", sep="/"))
    # pp_genes <- as.character(read.table(proto_file, header = FALSE,
    #                                     stringsAsFactors = FALSE)$V1)

    for (i in seq_along(gala_ds_batches)) {
        sce_obj <- gala_list[[i]]
        sce_features <- rownames(sce_obj)
        mtcg_features <- c(grep("ATMG", sce_features), grep("ATCG", sce_features))
        sce_obj <- addPerCellQC(sce_obj,
                                subsets = list(Mito = mtcg_features))
        sce_qc <- quickPerCellQC(colData(sce_obj),
                                 sub.fields = "subsets_Mito_percent")
        sce_obj <- sce_obj[, !sce_qc$discard]
        gala_list[[i]] <- logNormCounts(sce_obj)
    }
    common_genes <- Reduce(intersect, lapply(gala_list, rownames))

    for (i in seq_along(gala_list)) {
        sce_obj <- gala_list[[i]]
        gala_list[[i]] <- sce_obj[common_genes, ]
    }
    print("Constructed Objects.")
    mbat_out <- do.call(multiBatchNorm, gala_list)
    model_gene_vars <- lapply(mbat_out, function(sce_obj) {
        modelGeneVar(sce_obj)
    })
    combined_gene_vars <- do.call(combineVar, model_gene_vars)
    chosen_hvgs <- getTopHVGs(combined_gene_vars, n = 5000)
    print("Generate HVGs.")
    combined_sce <- do.call(correctExperiments,
                            append(mbat_out, list(PARAM = NoCorrectParam())))
    combined_sce <- runPCA(combined_sce, subset_row=chosen_hvgs)
    combined_sce <- runUMAP(combined_sce, dimred="PCA")
    #combined_sce <- runTSNE(combined_sce, dimred="PCA")
    print("Completed correctExperiments")
    fastmnn_out <- do.call(fastMNN,
                           append(mbat_out, list(subset.row = chosen_hvgs)))
    fastmnn_out <- runUMAP(fastmnn_out, dimred = "corrected")

    print("Completed FastMNN.")
    fmnn_adata <- zellkonverter::SCE2AnnData(fastmnn_out)
    hd5out_file <- paste(out_dir, "fastmnn_gala.h5ad", sep = "/")
    anndata::write_h5ad(fmnn_adata, filename = hd5out_file)

    rdsout_file <- paste(out_dir, "fastmnn_gala.rds", sep = "/")
    saveRDS(fastmnn_out, rdsout_file)
    print("Saved FastMNN objects.")
}

run_jb_gala_flt2_dset <- function(data_dir, proto_file, out_dir) {
    # tic("PREP DATASET TIME")
    start.time <- Sys.time()
    gala_file <- paste(data_dir, "GSE158761/GSE158761-FLT1.h5ad", sep = "/")
    gala_ds <- read_h5ad(gala_file)
    gala_ds_matrix <- gala_ds$X
    gala_ds_obs <- gala_ds$obs
    gala_ds_counts <- transpose_dgRMatrix(gala_ds_matrix)
    gala_ds_batches <- levels(gala_ds_obs[, "batch"])
    gala_ds_list <- lapply(gala_ds_batches, function(bx) {
        gala_ds_bcodes <- gala_ds_obs[gala_ds_obs["batch"] == bx, "Assay"]
        gala_ds_bcounts <-  gala_ds_counts[, gala_ds_bcodes]
        gala_ds_bcounts %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            as.SingleCellExperiment()
    })
    names(gala_ds_list) <- gala_ds_batches
    jb_file <- paste(data_dir,  "E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT2.h5ad",
                     sep = "/")
    jb_ds <- read_h5ad(jb_file)
    jb_ds_matrix <- jb_ds$X
    jb_ds_obs <- jb_ds$obs
    jb_ds_counts <- transpose_dgRMatrix(jb_ds_matrix)
    jb_ds_batches <- levels(jb_ds_obs[, "batch"])
    jb_ds_list <- lapply(jb_ds_batches, function(bx) {
        jb_ds_bcodes <- jb_ds_obs[jb_ds_obs["batch"] == bx, "Assay"]
        jb_ds_bcounts <-  jb_ds_counts[, jb_ds_bcodes]
        jb_ds_bcounts %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            as.SingleCellExperiment()
    })
    names(jb_ds_list) <- jb_ds_batches
    jb_gala_list <- c(gala_ds_list, jb_ds_list)
    jb_gala_batches <- c(gala_ds_batches, jb_ds_batches)
    saveRDS(jb_gala_list, file = paste(out_dir, "fastmnn_jb_gala_flt2_list.rds",
                                       sep = "/"))
    # pp_genes <- as.character(read.table(proto_file, header = FALSE,
    #                                     stringsAsFactors = FALSE)$V1)
    for (i in seq_along(jb_gala_batches)) {
        sce_obj <- jb_gala_list[[i]]
        sce_features <- rownames(sce_obj)
        mtcg_features <- c(grep("ATMG", sce_features), grep("ATCG", sce_features))
        sce_obj <- addPerCellQC(sce_obj,
                                subsets = list(Mito = mtcg_features))
        sce_qc <- quickPerCellQC(colData(sce_obj),
                                 sub.fields = "subsets_Mito_percent")
        sce_obj <- sce_obj[, !sce_qc$discard]
        #sce_obj <- sce_obj[, -mtcg_features]
        jb_gala_list[[i]] <- logNormCounts(sce_obj)
    }
    common_genes <- Reduce(intersect, lapply(jb_gala_list, rownames))

    for (i in seq_along(jb_gala_list)) {
        sce_obj <- jb_gala_list[[i]]
        jb_gala_list[[i]] <- sce_obj[common_genes, ]
        print(jb_gala_list[[i]])
    }
    print("Constructed Objects.")
    #toc()
    #tic("INTEGRATION TIME")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("PREP TIME : ", time.taken, units(time.taken),  "\n")
    start.time <- Sys.time()
    mbat_out <- do.call(multiBatchNorm, jb_gala_list)
    model_gene_vars <- lapply(mbat_out, function(sce_obj) {
        modelGeneVar(sce_obj)
    })
    combined_gene_vars <- do.call(combineVar, model_gene_vars)
    chosen_hvgs <- getTopHVGs(combined_gene_vars, var.field="p.value", n=25000)
    print("Generated HVGs.")
    combined_sce <- do.call(correctExperiments,
                            append(mbat_out, list(PARAM = NoCorrectParam())))
    #combined_sce <- runPCA(combined_sce, subset_row = chosen_hvgs)
    #combined_sce <- runUMAP(combined_sce, dimred = "PCA")
    #combined_sce <- runTSNE(combined_sce, dimred="PCA")
    #print("Completed correctExperiments")
    fastmnn_out <- do.call(fastMNN,
                           append(mbat_out, list(subset.row = chosen_hvgs)))
    fastmnn_out <- runUMAP(fastmnn_out, dimred = "corrected")
    print("Completed FastMNN.")
    #toc()
    #tic("POST-INTEGRATION TIME")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("INTEGRATION TIME : ", time.taken, units(time.taken),  "\n")
    start.time <- Sys.time()
    h5srt_file <- paste(out_dir, "fastmnn_jb_gala_flt2.h5Seurat", sep = "/")
    # hd5out_file <- paste(out_dir, "fastmnn_jb_gala_flt2.h5ad", sep = "/")
    reducedDims(combined_sce) <- reducedDims(fastmnn_out)
    reducedDimNames(combined_sce) <- c("corrected", "FUMAP")
    fmnn_srt <- as.Seurat(combined_sce)
    SaveH5Seurat(fmnn_srt, h5srt_file, overwrite = TRUE)
    Convert(h5srt_file, "h5ad", overwrite = TRUE)
    #fmnn_adata <- zellkonverter::SCE2AnnData(fastmnn_out)
    #anndata::write_h5ad(fmnn_adata, filename = hd5out_file)
    #
    rdsout_file <- paste(out_dir, "fastmnn_jb_gala_flt2.rds", sep = "/")
    saveRDS(fastmnn_out, rdsout_file)
    print("Saved FastMNN objects.")
    #toc()
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("POST-INTEGRATION TIME : ", time.taken, units(time.taken), "\n")
}

run_jb_gala_dset_integration <- function(data_dir, proto_file, out_dir) {
    start.time <- Sys.time()
    gala_file <- paste(data_dir, "GSE158761/GSE158761.h5ad", sep = "/")
    gala_ds <- read_h5ad(gala_file)
    gala_ds_matrix <- gala_ds$X
    gala_ds_obs <- gala_ds$obs
    gala_ds_counts <- transpose_dgRMatrix(gala_ds_matrix)
    gala_ds_batches <- levels(gala_ds_obs[, "batch"])
    gala_ds_list <- lapply(gala_ds_batches, function(bx) {
        gala_ds_bcodes <- gala_ds_obs[gala_ds_obs["batch"] == bx, "Assay"]
        gala_ds_bcounts <-  gala_ds_counts[, gala_ds_bcodes]
        gala_ds_bcounts %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            as.SingleCellExperiment()
    })
    names(gala_ds_list) <- gala_ds_batches
    jb_file <- paste(data_dir,  "E-GEOD-121619/E-GEOD-121619-FULL-LB.h5ad",
                     sep = "/")
    jb_ds <- read_h5ad(jb_file)
    jb_ds_matrix <- jb_ds$X
    jb_ds_obs <- jb_ds$obs
    jb_ds_counts <- transpose_dgRMatrix(jb_ds_matrix)
    jb_ds_batches <- levels(jb_ds_obs[, "batch"])
    jb_ds_list <- lapply(jb_ds_batches, function(bx) {
        jb_ds_bcodes <- jb_ds_obs[jb_ds_obs["batch"] == bx, "Assay"]
        jb_ds_bcounts <-  jb_ds_counts[, jb_ds_bcodes]
        jb_ds_bcounts %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            as.SingleCellExperiment()
    })
    names(jb_ds_list) <- jb_ds_batches
    jb_gala_list <- c(gala_ds_list, jb_ds_list)
    jb_gala_batches <- c(gala_ds_batches, jb_ds_batches)
    saveRDS(jb_gala_list, file = paste(out_dir, "fastmnn_jb_gala_list.rds",
                                       sep = "/"))
    # pp_genes <- as.character(read.table(proto_file, header = FALSE,
    #                                     stringsAsFactors = FALSE)$V1)
    for (i in seq_along(jb_gala_batches)) {
        sce_obj <- jb_gala_list[[i]]
        sce_features <- rownames(sce_obj)
        mtcg_features <- c(grep("ATMG", sce_features), grep("ATCG", sce_features))
        sce_obj <- addPerCellQC(sce_obj,
                                subsets = list(Mito = mtcg_features))
        sce_qc <- quickPerCellQC(colData(sce_obj),
                                 sub.fields = "subsets_Mito_percent")
        sce_obj <- sce_obj[, !sce_qc$discard]
        # sce_obj <- sce_obj[, -mtcg_features]
        jb_gala_list[[i]] <- logNormCounts(sce_obj)
    }
    common_genes <- Reduce(intersect, lapply(jb_gala_list, rownames))

    for (i in seq_along(jb_gala_list)) {
        sce_obj <- jb_gala_list[[i]]
        jb_gala_list[[i]] <- sce_obj[common_genes, ]
        print(jb_gala_list[[i]])
    }
    print("Constructed Objects.")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("PREP TIME : ", time.taken, units(time.taken),  "\n")
    start.time <- Sys.time()
    mbat_out <- do.call(multiBatchNorm, jb_gala_list)
    model_gene_vars <- lapply(mbat_out, function(sce_obj) {
        modelGeneVar(sce_obj)
    })
    combined_gene_vars <- do.call(combineVar, model_gene_vars)
    chosen_hvgs <- getTopHVGs(combined_gene_vars, var.field="p.value", n=5000)
    print("Generated HVGs.")
    combined_sce <- do.call(correctExperiments,
                            append(mbat_out, list(PARAM = NoCorrectParam())))
    #combined_sce <- runPCA(combined_sce, subset_row = chosen_hvgs)
    #combined_sce <- runUMAP(combined_sce, dimred = "PCA")
    #combined_sce <- runTSNE(combined_sce, dimred="PCA")
    #print("Completed correctExperiments")
    fastmnn_out <- do.call(fastMNN,
                           append(mbat_out, list(subset.row = chosen_hvgs)))
    fastmnn_out <- runUMAP(fastmnn_out, dimred = "corrected")
    print("Completed FastMNN.")
    #toc()
    #tic("POST-INTEGRATION TIME")
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("INTEGRATION TIME : ", time.taken, units(time.taken),  "\n")
    start.time <- Sys.time()
    h5srt_file <- paste(out_dir, "fastmnn_jb_gala.h5Seurat", sep = "/")
    hd5out_file <- paste(out_dir, "fastmnn_jb_gala.h5ad", sep = "/")
    reducedDims(combined_sce) <- reducedDims(fastmnn_out)
    reducedDimNames(combined_sce) <- c("corrected", "FUMAP")
    fmnn_srt <- as.Seurat(combined_sce)
    SaveH5Seurat(fmnn_srt, h5srt_file, overwrite = TRUE)
    Convert(h5srt_file, "h5ad", overwrite = TRUE)
    #fmnn_adata <- zellkonverter::SCE2AnnData(fastmnn_out)
    #anndata::write_h5ad(fmnn_adata, filename = hd5out_file)
    #
    rdsout_file <- paste(out_dir, "fastmnn_jb_gala.rds", sep = "/")
    saveRDS(fastmnn_out, rdsout_file)
    print("Saved FastMNN objects.")
    #toc()
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("POST-INTEGRATION TIME : ", time.taken, units(time.taken), "\n")
}


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) == 0){
    cat("Running default \n")
    run_jb_gala_flt2_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
} else {
    if (args[1] == "gala") {
        run_gala_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    
    if (args[1] == "jb_gala_flt2") {
        run_jb_gala_flt2_dset(DATA_DIR, PROTPLAST_GENE_LIST,
                              ANAYSIS_OUT_DIR)
    }
    
    if (args[1] == "jb_gala") {
        run_jb_gala_dset_integration(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    
    if (args[1] == "default") {
        run_jb_gala_flt2_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
}
