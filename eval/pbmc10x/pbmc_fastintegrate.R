suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(Seurat))
suppressMessages(library(plotly))
suppressMessages(library(grid))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(pbmcapply))
suppressMessages(library(stringr))
suppressMessages(library(FastIntegration))

args <- commandArgs(trailingOnly=TRUE)
cat("Running with args : ", args, "\n")
ndata <- as.integer(args[1])
cat("Setting No. of datasets : ", ndata, "\n")
tmp_base_dir <- "../../data/FastIntegrate"
seurat_obj_dir <- "../../data/FastIntegrate/Seurat_Objects/"

if (ndata == 2) {
    use.sample <- c("n06", "n07")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate2/")
}

if (ndata == 3) {
    use.sample <- c("n02", "n06", "n07")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate3/")
}

if (ndata == 5) {
    use.sample <- c("n01", "n02", "n06", "n07", "n08")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate5/")
}

if (ndata == 9) {
    use.sample <- c("n01", "n02", "n06", "n07", "n08",
                    "n03", "n05", "n13", "n14")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate9/")
}

if (ndata == 13) {
    use.sample <- c("n01", "n02", "n06", "n07", "n08", "n03",
                    "n05", "n09", "n10", "n11", "n12", "n13", "n14")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate13/")
}

if (ndata == 16) {
    use.sample <- c("n01", "n02", "n06", "n07", "n08", "n03", "n05", "n09",
                    "n10", "n11", "n12", "n13", "n14", "n70", "n71", "n74")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate16/")
}

if (ndata == 17) {
    use.sample <- c("n01", "n02", "n06", "n07", "n08", "n03", "n05", "n09", "n10",
                    "n11", "n12", "n13", "n14", "n70", "n71", "n74", "n15")
    temp.dir <- paste0(tmp_base_dir, "/FastIntegrate17/")
}

ncores <- 2

cat("Setting sample names : ", use.sample, "\n")

read_seu <- function(dir, sample.name) {
    seu <- readRDS(dir)
    return(seu)
}

list.filenames <- list.files(path=seurat_obj_dir, pattern=".rds$") %>%
    .[match(use.sample, gsub("_seurat.rds", "", .))]

rna.list <- list()

for (i in 1:length(list.filenames)) {
    rna.list[[i]]<-read_seu(dir = paste0(seurat_obj_dir, list.filenames[i]),
                            sample.name = use.sample[i])
}

#names(rna.list) <- list.filenames %>% gsub("_seurat.rds","",.)

prep_start_time <- Sys.time()
overlapped.gene <- Reduce(intersect, lapply(rna.list, rownames))
for (i in 1:length(rna.list)) {
    rna.list[[i]] <- subset(rna.list[[i]], features = overlapped.gene)
    rna.list[[i]] <- NormalizeData(rna.list[[i]])
    rna.list[[i]] <- FindVariableFeatures(rna.list[[i]])
    rna.list[[i]] <- RenameCells(rna.list[[i]],
                                 new.names = paste0(Cells(rna.list[[i]]), "--", i))
}
prep_end_time <- Sys.time()
prep_runtime <- prep_end_time - prep_start_time
cat("Overlapped genes size : ", length(overlapped.gene), "\n")
cat("Subset Features and Normalize Data Runtime is ", prep_runtime,
    units(prep_runtime), "\n")

prep_start_time <- Sys.time()
BuildIntegrationFile(rna.list = rna.list, tmp.dir = temp.dir, nCores = ncores)
FastFindAnchors(tmp.dir = temp.dir, nCores = ncores)
genes <- readRDS(paste0(temp.dir, "/FastIntegrationTmp/raw/1.rds"))
genes <- rownames(genes)
idx <- split(1:length(genes), cut(1:length(genes), ncores, labels = FALSE))
prep_end_time <- Sys.time()
prep_runtime <- prep_end_time - prep_start_time
cat("Gene size begin split : ", length(genes), "\n")
cat("BuildIntegrationFile and FastFindAnchors Runtime is ", prep_runtime,
    units(prep_runtime), "\n")

prep_start_time <- Sys.time()
pbmclapply(
    1:ncores, function(i) {
        rna.integrated <- FastIntegration(tmp.dir = temp.dir, npcs = 1:30,
                                          slot = "data",
                                          features.to.integrate = genes[idx[[i]]])
        saveRDS(rna.integrated,
                paste0(temp.dir, "/FastIntegrationTmp/inte/inte_pbmc", i, ".rds"),
                compress=FALSE)
    }, mc.cores = ncores
)
prep_end_time <- Sys.time()
prep_runtime <- prep_end_time - prep_start_time
cat("FastIntegration Runtime is ", prep_runtime, units(prep_runtime), "\n")
print(rna.integrated)
