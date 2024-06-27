library(dplyr)
library(anndata)
library(Seurat)
library(SeuratDisk)
library(tictoc)

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
    tic("PREP DATASET TIME")
    # start.time <- Sys.time()
    gala_file <- paste(data_dir, "GSE158761/GSE158761.h5ad", sep = "/")
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
            SCTransform(variable.features.n = 20000)
    })
    saveRDS(gala_list, file = paste(out_dir, "seurat_gala_list.rds", sep = "/"))
    toc()
    tic("INTEGRATION TIME")
    pp_genes <- as.character(read.table(proto_file, header = FALSE,
                                        stringsAsFactors = FALSE)$V1)
    gala_features <- SelectIntegrationFeatures(object.list = gala_list,
                                               nfeatures = 25000)
    gala_features <- gala_features[-c(grep("ATMG", gala_features),
                                      grep("ATCG", gala_features),
                                      sort(match(pp_genes, gala_features)))]
    gala_list <- PrepSCTIntegration(object.list = gala_list,
                                    anchor.features = gala_features,
                                    verbose = TRUE)
    gala_anchors <- FindIntegrationAnchors(object.list = gala_list,
                                           normalization.method = "SCT",
                                           anchor.features = gala_features,
                                           verbose = TRUE)
    gala_integrated <- IntegrateData(anchorset = gala_anchors,
                                     normalization.method = "SCT",
                                     verbose = TRUE)
    gala_integrated <- RunPCA(gala_integrated, verbose = TRUE)
    gala_integrated <- RunUMAP(gala_integrated, reduction = "pca", dims = 1:30)
    saveRDS(gala_integrated, paste(out_dir, "seurat_integrated_ath_gala.rds",
                                   sep = "/"))
    toc()
    tic("POSTPP TIME")
    gala_integ_matrix <- gala_integrated@assays$integrated@data
    gala_integ_pca <- gala_integrated@reductions$pca@cell.embeddings
    gala_integ_umap <- gala_integrated@reductions$umap@cell.embeddings

    gala_data_file <- paste(out_dir, "scanpy_normalized_ath_gala.h5ad",
                            sep = "/")
    gala_ad_object <- read_h5ad(gala_data_file)
    gala_meta_data <- gala_ad_object$obs
    gala_seurat_results <- list(
        "integ_norm" = gala_integ_matrix,
        "meta_data" = gala_meta_data,
        "seurat_embeddings" = gala_integ_pca,
        "seurat_umap" = gala_integ_umap
    )
    saveRDS(gala_seurat_results,
            paste(out_dir, "seurat_ath_gala.rds", sep = "/"))
    SaveH5Seurat(gala_integrated,
                 filename = paste(out_dir, "seurat_ath_gala.h5Seurat",
                                  sep = "/"),
                 overwrite = TRUE)
    Convert(paste(out_dir, "seurat_ath_gala.h5Seurat", sep = "/"),
            "h5ad", overwrite = TRUE)
    toc()
}

run_jb_gala_flt1_dset <- function(data_dir, proto_file, out_dir) {
    tic("PREP DATASET TIME")
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
            SCTransform(variable.features.n = 20000)
    })
    names(gala_ds_list) <- gala_ds_batches
    jb_file <- paste(data_dir,  "E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT1.h5ad",
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
            SCTransform(variable.features.n = 20000)
    })
    names(jb_ds_list) <- jb_ds_batches
    jb_gala_list <- c(gala_ds_list, jb_ds_list)
    saveRDS(jb_gala_list, file = paste(out_dir, "seurat_jb_gala_flt1_list.rds",
                                       sep = "/"))
    toc()
    tic("INTEGRATION TIME")
    pp_genes <- as.character(read.table(proto_file, header = FALSE,
                                        stringsAsFactors = FALSE)$V1)
    jb_gala_features <- SelectIntegrationFeatures(object.list = jb_gala_list,
                                                  nfeatures = 25000)
    jb_gala_features <- jb_gala_features[-c(grep("ATMG", jb_gala_features),
                                            grep("ATCG", jb_gala_features),
                                            sort(match(pp_genes,
                                                       jb_gala_features)))]
    jb_gala_list <- PrepSCTIntegration(object.list = jb_gala_list,
                                       anchor.features = jb_gala_features,
                                       verbose = TRUE)
    jb_gala_anchors <- FindIntegrationAnchors(object.list = jb_gala_list,
                                              normalization.method = "SCT",
                                              anchor.features = jb_gala_features,
                                              verbose = TRUE)
    jb_gala_integrated <- IntegrateData(anchorset = jb_gala_anchors,
                                        normalization.method = "SCT",
                                        verbose = TRUE)
    jb_gala_integrated <- RunPCA(jb_gala_integrated, verbose = TRUE)
    jb_gala_integrated <- RunUMAP(jb_gala_integrated, reduction = "pca",
                                  dims = 1:30)
    saveRDS(jb_gala_integrated,
            paste(out_dir, "seurat_integrated_ath_jb_gala_flt1.rds", sep = "/"))
    toc()
    tic("POSTPP TIME")
    jb_gala_integ_matrix <- jb_gala_integrated@assays$integrated@data
    jb_gala_integ_pca <- jb_gala_integrated@reductions$pca@cell.embeddings
    jb_gala_integ_umap <- jb_gala_integrated@reductions$umap@cell.embeddings

    jb_gala_data_file <- paste(out_dir,
                               "scanpy_normalized_ath_jb_gala_flt1.h5ad",
                               sep = "/")
    jb_gala_ad_object <- read_h5ad(jb_gala_data_file)
    jb_gala_meta_data <- jb_gala_ad_object$obs
    jb_gala_seurat_results <- list(
        "integ_norm" = jb_gala_integ_matrix,
        "meta_data" = jb_gala_meta_data,
        "seurat_embeddings" = jb_gala_integ_pca,
        "seurat_umap" = jb_gala_integ_umap
    )
    saveRDS(jb_gala_seurat_results,
            paste(out_dir, "seurat_ath_jb_gala_flt1.rds", sep = "/"))
    SaveH5Seurat(jb_gala_integrated,
                 filename = paste(out_dir, "seurat_ath_jb_gala_flt1.h5Seurat",
                                  sep = "/"),
                 overwrite = TRUE)
    Convert(paste(out_dir, "seurat_ath_jb_gala_flt1.h5Seurat", sep = "/"),
            "h5ad", overwrite = TRUE)
    toc()
}


run_jb_gala_flt2_dset <- function(data_dir, proto_file, out_dir) {
    tic("PREP DATASET TIME")
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
            SCTransform(variable.features.n = 20000)
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
            SCTransform(variable.features.n = 20000)
    })
    names(jb_ds_list) <- jb_ds_batches
    jb_gala_list <- c(gala_ds_list, jb_ds_list)
    saveRDS(jb_gala_list, file = paste(out_dir, "seurat_jb_gala_flt2_list.rds",
                                       sep = "/"))
    toc()
    tic("INTEGRATION TIME")
    pp_genes <- as.character(read.table(proto_file, header = FALSE,
                                        stringsAsFactors = FALSE)$V1)
    jb_gala_features <- SelectIntegrationFeatures(object.list = jb_gala_list,
                                                  nfeatures = 25000)
    jb_gala_features <- jb_gala_features[-c(grep("ATMG", jb_gala_features),
                                            grep("ATCG", jb_gala_features),
                                            sort(match(pp_genes,
                                                       jb_gala_features)))]
    jb_gala_list <- PrepSCTIntegration(object.list = jb_gala_list,
                                       anchor.features = jb_gala_features,
                                       verbose = TRUE)
    jb_gala_anchors <- FindIntegrationAnchors(object.list = jb_gala_list,
                                              normalization.method = "SCT",
                                              anchor.features = jb_gala_features,
                                              verbose = TRUE)
    jb_gala_integrated <- IntegrateData(anchorset = jb_gala_anchors,
                                        normalization.method = "SCT",
                                        verbose = TRUE)
    jb_gala_integrated <- RunPCA(jb_gala_integrated, verbose = TRUE)
    jb_gala_integrated <- RunUMAP(jb_gala_integrated, reduction = "pca",
                                  dims = 1:30)
    saveRDS(jb_gala_integrated,
            paste(out_dir, "seurat_integrated_ath_jb_gala_flt2.rds", sep = "/"))
    toc()
    tic("POSTPP TIME")
    jb_gala_integ_matrix <- jb_gala_integrated@assays$integrated@data
    jb_gala_integ_pca <- jb_gala_integrated@reductions$pca@cell.embeddings
    jb_gala_integ_umap <- jb_gala_integrated@reductions$umap@cell.embeddings

    jb_gala_data_file <- paste(out_dir,
                               "scanpy_normalized_ath_jb_gala_flt2.h5ad",
                               sep = "/")
    jb_gala_ad_object <- read_h5ad(jb_gala_data_file)
    jb_gala_meta_data <- jb_gala_ad_object$obs
    jb_gala_seurat_results <- list(
        "integ_norm" = jb_gala_integ_matrix,
        "meta_data" = jb_gala_meta_data,
        "seurat_embeddings" = jb_gala_integ_pca,
        "seurat_umap" = jb_gala_integ_umap
    )
    saveRDS(jb_gala_seurat_results,
            paste(out_dir, "seurat_ath_jb_gala_flt2.rds", sep = "/"))
    SaveH5Seurat(jb_gala_integrated,
                 filename = paste(out_dir, "seurat_ath_jb_gala_flt2.h5Seurat",
                                  sep = "/"),
                 overwrite = TRUE)
    Convert(paste(out_dir, "seurat_ath_jb_gala_flt2.h5Seurat", sep = "/"),
            "h5ad", overwrite = TRUE)
    toc()
}




run_jb_gala_dset <- function(data_dir, proto_file, out_dir) {
    tic("PREP DATASET TIME")
    gala_file <- paste(data_dir, "GSE158761/GSE158761.h5ad", sep = "/")
    gala_ds <- read_h5ad(gala_file)
    gala_ds_matrix <- gala_ds$X
    gala_ds_obs <- gala_ds$obs
    gala_ds_counts <- transpose_dgRMatrix(gala_ds_matrix)
    gala_ds_batches <- levels(gala_ds_obs[, "batch"])
    gala_ds_list <- lapply(gala_ds_batches, function(bx) {
        gala_ds_bcodes <- gala_ds_obs[gala_ds_obs["batch"] == bx,
                                      "Assay"]
        gala_ds_bcounts <-  gala_ds_counts[, gala_ds_bcodes]
        gala_ds_bcounts %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            SCTransform(variable.features.n = 20000)
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
            SCTransform(variable.features.n = 20000)
    })
    names(jb_ds_list) <- jb_ds_batches
    jb_gala_list <- c(gala_ds_list, jb_ds_list)
    saveRDS(jb_gala_list, file = paste(out_dir, "seurat_jb_gala_list.rds",
                                       sep = "/"))
    toc()
    tic("INTEGRATION TIME")
    pp_genes <- as.character(read.table(proto_file, header = FALSE,
                                        stringsAsFactors = FALSE)$V1)

    jb_gala_features <- SelectIntegrationFeatures(object.list = jb_gala_list,
                                                  nfeatures = 25000)
    jb_gala_features <- jb_gala_features[-c(grep("ATMG", jb_gala_features),
                                            grep("ATCG", jb_gala_features),
                                            sort(match(pp_genes,
                                                       jb_gala_features)))]
    jb_gala_list <- PrepSCTIntegration(object.list = jb_gala_list,
                                       anchor.features = jb_gala_features,
                                       verbose = TRUE)
    jb_gala_anchors <- FindIntegrationAnchors(object.list = jb_gala_list,
                                              normalization.method = "SCT",
                                              anchor.features = jb_gala_features,
                                              verbose = TRUE)
    jb_gala_integrated <- IntegrateData(anchorset = jb_gala_anchors,
                                        normalization.method = "SCT",
                                        verbose = TRUE)
    jb_gala_integrated <- RunPCA(jb_gala_integrated, verbose = TRUE)
    jb_gala_integrated <- RunUMAP(jb_gala_integrated, reduction = "pca",
                                  dims = 1:30)
    saveRDS(jb_gala_integrated,
            paste(out_dir, "seurat_integrated_ath_jb_gala.rds", sep = "/"))
    toc()
    tic("POSTPP TIME")
    jb_gala_integ_matrix <- jb_gala_integrated@assays$integrated@data
    jb_gala_integ_pca <- jb_gala_integrated@reductions$pca@cell.embeddings
    jb_gala_integ_umap <- jb_gala_integrated@reductions$umap@cell.embeddings

    jb_gala_data_file <- paste(out_dir, "scanpy_normalized_ath_jb_gala.h5ad",
                               sep = "/")
    jb_gala_ad_object <- read_h5ad(jb_gala_data_file)
    jb_gala_meta_data <- jb_gala_ad_object$obs
    jb_gala_seurat_results <- list(
        "integ_norm" = jb_gala_integ_matrix,
        "meta_data" = jb_gala_meta_data,
        "seurat_embeddings" = jb_gala_integ_pca,
        "seurat_umap" = jb_gala_integ_umap
    )
    saveRDS(jb_gala_seurat_results,
            paste(out_dir, "seurat_ath_jb_gala.rds", sep = "/"))
    SaveH5Seurat(jb_gala_integrated,
                 filename = paste(out_dir, "seurat_ath_jb_gala.h5Seurat",
                                  sep = "/"),
                 overwrite = TRUE)
    Convert(paste(out_dir, "seurat_ath_jb_gala.h5Seurat", sep = "/"),
            "h5ad", overwrite = TRUE)
    toc()
}


run_copilot_dset <- function(data_dir, proto_file, out_dir) {
    tic("PREP DATASET TIME")
    batches <- c("sc_1", "sc_31", "tnw2", "sc_11", "sc_51", "sc_37",
                 "sc_9_at", "tnw1", "sc_40", "sc_12", "col0", "sc_30",
                 "sc_10_at")
    cpt_rc_ads <- lapply(batches, function(bx) {
        h5ad_file <- paste(data_dir,  "GSE152766",
                           paste(bx, ".h5ad", sep = ""), sep = "/")
        read_h5ad(h5ad_file)
    })
    names(cpt_rc_ads) <- batches
    cpt_rc_list <- lapply(cpt_rc_ads, function(ad_object) {
        transpose_dgRMatrix(ad_object$X) %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            SCTransform(variable.features.n = 2000)
    })
    names(cpt_rc_list) <- batches
    toc()
    tic("INTEGRATION TIME")
    pp_genes <- as.character(read.table(proto_file, header = FALSE,
                                        stringsAsFactors = FALSE)$V1)
    cpt_rc_features <- SelectIntegrationFeatures(object.list = cpt_rc_list,
                                                 nfeatures = 2000)
    cpt_rc_features <- cpt_rc_features[-c(grep("ATMG", cpt_rc_features),
                                          grep("ATCG", cpt_rc_features),
                                          sort(match(pp_genes,
                                                     cpt_rc_features)))]
    cpt_rc_list <- PrepSCTIntegration(object.list = cpt_rc_list,
                                      anchor.features = cpt_rc_features,
                                      verbose = TRUE)
    cpt_rc_anchors <- FindIntegrationAnchors(object.list = cpt_rc_list,
                                             normalization.method = "SCT",
                                             anchor.features = cpt_rc_features,
                                             verbose = TRUE)
    cpt_rc_integrated <- IntegrateData(anchorset = cpt_rc_anchors,
                                       normalization.method = "SCT",
                                       verbose = TRUE)
    cpt_rc_integrated <- RunPCA(cpt_rc_integrated, verbose = TRUE)
    cpt_rc_integrated <- RunUMAP(cpt_rc_integrated, reduction = "pca",
                                 dims = 1:30)
    saveRDS(cpt_rc_integrated,
            paste(out_dir, "seurat_integrated_ath_copilot.rds", sep = "/"))
    toc()
    tic("POSTPP TIME")
    cpt_integ_matrix <- cpt_rc_integrated@assays$integrated@data
    cpt_integ_pca <- cpt_rc_integrated@reductions$pca@cell.embeddings
    cpt_integ_umap <- cpt_rc_integrated@reductions$umap@cell.embeddings

    cpt_data_file <- paste(out_dir, "scanpy_normalized_ath_copilot.h5ad",
                           sep = "/")
    cpt_ad_object <- read_h5ad(cpt_data_file)
    cpt_meta_data <- cpt_ad_object$obs
    cpt_seurat_results <- list(
        "integ_norm" = cpt_integ_matrix,
        "meta_data" = cpt_meta_data,
        "seurat_embeddings" = cpt_integ_pca,
        "seurat_umap" = cpt_integ_umap
    )
    saveRDS(cpt_seurat_results,
            paste(out_dir, "seurat_ath_copilot.rds", sep = "/"))

    SaveH5Seurat(cpt_rc_integrated,
                 filename = paste(out_dir, "seurat_ath_copilot.h5Seurat",
                                  sep = "/"),
                 overwrite = TRUE)
    Convert(paste(out_dir, "seurat_ath_copilot.h5Seurat", sep = "/"),
            "h5ad", overwrite = TRUE)
    toc()
}


run_wildtype <- function(data_dir, proto_file, out_dir) {
    tic("PREP DATASET TIME")
    t19_wt_ad <- read_h5ad(paste(data_dir,
                                 "E-GEOD-121619/E-GEOD-121619-LB.h5ad",
                                 sep = "/"))
    t13_wt_ad <- read_h5ad(paste(data_dir,
                                 "E-GEOD-123013/E-GEOD-123013-LB.h5ad",
                                 sep = "/"))
    t13_rhd_ad <- read_h5ad(paste(data_dir,
                                  "E-GEOD-123013/E-GEOD-123013-RHD.h5ad",
                                  sep = "/"))
    t61_wt_ad <- read_h5ad(paste(data_dir,
                                 "E-GEOD-158761/E-GEOD-158761-LB.h5ad",
                                 sep = "/"))
    t19_wt_counts <- transpose_dgRMatrix(t19_wt_ad$X)
    t13_wt_counts <- transpose_dgRMatrix(t13_wt_ad$X)
    t13_rhd_counts <- transpose_dgRMatrix(t13_rhd_ad$X)
    t61_wt_counts <- transpose_dgRMatrix(t61_wt_ad$X)
    #
    t19_wt_sobj <- t19_wt_counts %>%
        CreateSeuratObject(min.cells = 3, min.features = 200) %>%
        SCTransform(variable.features.n = 20000)
    saveRDS(t19_wt_sobj,
            paste(data_dir, "E-GEOD-121619/E-GEOD-121619-LB-SCT.rds",
                  sep = "/"))
    #
    t13_wt_sobj <- t13_wt_counts %>%
        CreateSeuratObject(min.cells = 3, min.features = 200) %>%
        SCTransform(variable.features.n = 20000)
    saveRDS(t13_wt_sobj,
            paste(data_dir, "E-GEOD-123013/E-GEOD-123013-LB-SCT.rds",
                  sep = "/"))
    #
    t13_rhd_sobj <- t13_rhd_counts %>%
        CreateSeuratObject(min.cells = 3, min.features = 200) %>%
        SCTransform(variable.features.n = 20000)
    saveRDS(t13_rhd_sobj,
            paste(data_dir, "E-GEOD-123013/E-GEOD-123013-RHD-SCT.rds",
                  sep = "/"))
    #
    t61_wt_sobj <- t61_wt_counts %>%
        CreateSeuratObject(min.cells = 3, min.features = 200) %>%
        SCTransform(variable.features.n = 20000)
    saveRDS(t61_wt_sobj,
            paste(data_dir, "E-GEOD-158761/E-GEOD-158761-LB-SCT.rds",
                  sep = "/"))
    #
    toc()
    tic("INTEGRATION TIME")
    pp_genes <- as.character(read.table(proto_file, header = FALSE,
                                        stringsAsFactors = FALSE)$V1)
    wt_rc_list <- list(
        "G121619" = t19_wt_sobj,
        "G123013" = t13_wt_sobj,
        "G158761"  = t61_wt_sobj
    )
    wt_rc_features <- SelectIntegrationFeatures(object.list = wt_rc_list,
                                                nfeatures = 25000)

    wt_rc_features <- wt_rc_features[-c(grep("ATMG", wt_rc_features),
                                        grep("ATCG", wt_rc_features),
                                        sort(match(pp_genes, wt_rc_features)))]

    wt_rc_list <- PrepSCTIntegration(object.list = wt_rc_list,
                                     anchor.features = wt_rc_features,
                                     verbose = TRUE)
    wt_rc_anchors <- FindIntegrationAnchors(object.list = wt_rc_list,
                                            normalization.method = "SCT",
                                            anchor.features = wt_rc_features,
                                            verbose = TRUE)
    wt_rc_integrated <- IntegrateData(anchorset = wt_rc_anchors,
                                      normalization.method = "SCT",
                                      verbose = TRUE)
    wt_rc_integrated <- RunPCA(wt_rc_integrated, verbose = TRUE)
    wt_rc_integrated <- RunUMAP(wt_rc_integrated, reduction = "pca",
                                dims = 1:30)
    saveRDS(wt_rc_integrated, paste(out_dir, "seurat_integrated_ath_wt.rds",
                                    sep = "/"))
    toc()
    tic("POSTPP TIME")
    wt_rc_combined <- merge(wt_rc_list[[1]],
                            list(wt_rc_list[[2]], wt_rc_list[[3]]),
                            project = "ATHWT")
    wt_merged_matrix <- wt_rc_combined@assays$SCT@data
    wt_integ_matrix <- wt_rc_integrated@assays$integrated@data
    wt_integ_pca <- wt_rc_integrated@reductions$pca@cell.embeddings
    wt_integ_umap <- wt_rc_integrated@reductions$umap@cell.embeddings
    wt_data_file <- paste(data_dir, "scanpy_ath_wt_raw.h5ad", sep = "/")
    wt_ad_object <- read_h5ad(wt_data_file)
    wt_exprs_raw <- as.matrix(wt_ad_object$T$X)
    wt_meta_data <- wt_ad_object$obs
    wt_seurat_results <- list(
        "exprs_norm" = wt_rc_combined,
        "integ_norm" = wt_integ_matrix,
        "meta_data" = wt_meta_data,
        "seurat_embeddings" = wt_integ_pca,
        "seurat_umap" = wt_integ_umap
    )
    saveRDS(wt_seurat_results, paste(out_dir, "seurat_ath_wt.rds", sep = "/"))

    wt_seurat_results_embs <- list(
        "meta_data" = wt_seurat_results$meta_data,
        "seurat_embeddings" = wt_seurat_results$seurat_embeddings,
        "seurat_umap" = wt_seurat_results$seurat_umap
    )
    saveRDS(wt_seurat_results_embs, paste(out_dir, "seurat_ath_wt_embeds.rds",
                                          sep = "/"))
    toc()
}


options(echo = TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) == 0){
    cat("Running default \n")
    run_jb_gala_flt2_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
} else {
    if (args[1] == "wt") {
        run_wildtype(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    if (args[1] == "copilot") {
        run_copilot_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    if (args[1] == "gala") {
        run_gala_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    if (args[1] == "jb_gala") {
        run_jb_gala_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    if (args[1] == "jb_gala_flt1") {
        run_jb_gala_flt1_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    if (args[1] == "jb_gala_flt2") {
        run_jb_gala_flt2_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
    if (args[1] == "default") {
        run_jb_gala_flt2_dset(DATA_DIR, PROTPLAST_GENE_LIST, ANAYSIS_OUT_DIR)
    }
}
