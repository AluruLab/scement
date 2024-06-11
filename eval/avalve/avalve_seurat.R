library(dplyr)
library(anndata)
library(Seurat)
library(SeuratDisk)
library(tictoc)

# Generate analyses with COMBAT
# source data is from DATA_DIR, and output is ANALYSIS_RESULTS_DIR
DATA_DIR <- "./aortic_valve/"
ANALYSIS_RESULTS_DIR <- "./avalve/"

transpose_dgRMatrix <- function(inmat) {
    if (class(inmat) != "dgRMatrix") {
        stop("inmat is not of class dgRMatrix")
    }
    out <- new("dgCMatrix",
        i = inmat@j, p = inmat@p,
        x = inmat@x, Dim = rev(inmat@Dim),
        Dimnames = rev(inmat@Dimnames)
    )
    out
}

integrate_seurat <- function(data_dir, out_dir) {
    tic("PREP DATASET TIME")
    batches <- c("H1", "H2", "C3", "C4")
    avalve_rc_ads <- lapply(batches, function(bx) {
        h5ad_file <- paste(data_dir,
            paste(bx, ".h5ad", sep = ""),
            sep = "/"
        )
        read_h5ad(h5ad_file)
    })
    names(avalve_rc_ads) <- batches
    avalve_rc_list <- lapply(avalve_rc_ads, function(ad_object) {
        transpose_dgRMatrix(ad_object$X) %>%
            CreateSeuratObject(min.cells = 3, min.features = 200) %>%
            SCTransform(variable.features.n = 20000)
    })
    names(avalve_rc_list) <- batches
    toc()
    tic("INTEGRATION TIME")
    avalve_rc_features <- SelectIntegrationFeatures(
        object.list = avalve_rc_list, nfeatures = 25000
    )
    #    avalve_rc_features <- avalve_rc_features[
    #                              -c(grep("ATMG", avalve_rc_features),
    #                               grep("ATCG", avalve_rc_features),
    #                               sort(match(pp_genes, avalve_rc_features)))]
    avalve_rc_list <- PrepSCTIntegration(
        object.list = avalve_rc_list,
        anchor.features = avalve_rc_features,
        verbose = TRUE
    )
    avalve_rc_anchors <- FindIntegrationAnchors(
        object.list = avalve_rc_list,
        normalization.method = "SCT",
        anchor.features = avalve_rc_features,
        verbose = TRUE
    )
    avalve_rc_integrated <- IntegrateData(
        anchorset = avalve_rc_anchors,
        normalization.method = "SCT",
        verbose = TRUE
    )
    avalve_rc_integrated <- RunPCA(avalve_rc_integrated, verbose = TRUE)
    avalve_rc_integrated <- RunUMAP(avalve_rc_integrated,
        reduction = "pca",
        dims = 1:30
    )
    saveRDS(
        avalve_rc_integrated,
        paste(out_dir, "seurat_integrated_avalve.rds", sep = "/")
    )
    toc()
    tic("POSTPP TIME")
    avalve_integ_matrix <- avalve_rc_integrated@assays$integrated@data
    avalve_integ_pca <- avalve_rc_integrated@reductions$pca@cell.embeddings
    avalve_integ_umap <- avalve_rc_integrated@reductions$umap@cell.embeddings

    avalve_data_file <- paste(out_dir, "scanpy_normalized_avalve.h5ad",
        sep = "/"
    )
    avalve_ad_object <- read_h5ad(avalve_data_file)
    avalve_meta_data <- avalve_ad_object$obs
    avalve_seurat_results <- list(
        "integ_norm" = avalve_integ_matrix,
        "meta_data" = avalve_meta_data,
        "seurat_embeddings" = avalve_integ_pca,
        "seurat_umap" = avalve_integ_umap
    )
    saveRDS(
        avalve_seurat_results,
        paste(out_dir, "seurat_avalve.rds", sep = "/")
    )

    SaveH5Seurat(avalve_rc_integrated,
        filename = paste(out_dir, "seurat_avalve.h5Seurat", sep = "/"),
        overwrite = TRUE
    )
    Convert(paste(out_dir, "seurat_avalve.h5Seurat", sep = "/"), "h5ad",
        overwrite = TRUE
    )
    toc()
}

options(echo = TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
if (args[1] == "seurat") {
    integrate_seurat(DATA_DIR, ANALYSIS_RESULTS_DIR)
}
