library(Seurat)
library(BPCells)
library(ggplot2)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# patients data
files_dir <- "../../data/covid_atlas/"
merged_dir <- "../../data/covid_atlas/"
meta_data_csv_file <- "../../data/covid_atlas/covid_metadata.csv"
cat("Merges data from samples of h5ad_full data", "\n")
cat("IN : ", files_dir, "OUT : ", merged_dir,  "\n")

patients_500k <- c(
  "S-S077", "S-M010-3", "S-M011", "S-S045", "S-M009-5", "S-S045",
  "S-M019", "S-M066", "S-S077", "S-M039-2", "S-M031-1", "S-M036-1",
  "S-S073-1", "S-M029", "S-M062-2", "S-M032", "S-M008-2", "S-M059-1",
  "S-S022-5", "S-M016", "S-S014", "S-M022", "S-M060-2", "S-M048",
  "S-S088-2", "S-M066", "S-M071", "S-S074-2", "S-S017", "S-M051",
  "S-M040-1", "S-M032", "S-M010-6", "S-M072", "S-M004-4", "S-S070-3",
  "S-S020", "S-S091", "S-M035-1", "S-M076-2", "S-S084", "S-S076-1",
  "S-M010-1", "S-S035-2", "S-S025", "S-M062-2", "S-S022-1", "S-M006",
  "S-S045", "S-M054", "S-S067", "S-M009-5", "S-S041", "S-S029",
  "S-M009-2", "S-S066", "S-S051", "S-S057", "S-M048", "S-M061-2",
  "S-M019", "S-S035-1", "S-S039", "S-M040-1", "S-M037", "S-S021-3",
  "S-M018", "S-M058-1", "S-S090-2", "S-S061", "S-M001", "S-M071",
  "S-S064", "S-M040-1", "S-M054", "S-M005", "S-S089-2", "S-S086-2",
  "S-M056", "S-M026-1", "S-S078", "S-S034", "S-M032", "S-M060-2",
  "S-S087-2"
)
patients_full <- c(
  "S-S070-2", "S-S070-3", "S-S069-3", "S-M056", "S-M044-1", "S-M043-1",
  "S-M048", "S-M044-2", "S-M043-2", "S-S054", "S-S056", "S-M042-1", "S-M041-1",
  "S-M049", "S-M046", "S-M047", "S-S055", "S-S057", "S-M045", "S-M041-2",
  "S-M042-2", "S-M055", "S-M053", "S-S067", "S-S065", "S-M051", "S-S064",
  "S-M054", "S-M052", "S-S068", "S-S066", "S-S059", "S-S060", "S-M050",
  "S-S061", "S-S062", "S-S063", "S-M061-1", "S-M061-2", "S-S073-1",
  "S-S074-1", "S-S074-2", "S-M062-2", "S-M058-1", "S-M058-2", "S-M063-1",
  "S-M063-2", "S-M059-1", "S-M059-2", "S-M060-1", "S-M060-2", "S-S075-1",
  "S-S076-1", "S-S076-2", "S-S035-1", "S-S035-2", "S-S035-3",
  "S-S035-4", "S-S036-1", "S-S036-2", "S-S036-3", "S-M064", "S-M066", "S-M067",
  "S-M068", "S-S091", "S-S090-2", "S-S092", "S-M076-2", "S-M074-2", "S-S088-2",
  "S-S085-2", "S-S089-2", "S-S087-2", "S-S086-2", "S-M077", "S-M078", "S-M079",
  "S-S024", "S-S026", "S-M013", "S-M014", "S-M015", "S-M016", "S-S023",
  "S-S025", "S-S027", "S-S034", "S-M025", "S-S033", "S-M028", "S-M029",
  "S-M027", "S-M026-1", "S-S032-3", "S-M026-2", "S-M026-3", "S-M004-1",
  "S-M004-2", "S-M004-3", "S-S022-1", "S-S022-2", "S-S022-3", "S-S022-4",
  "S-S022-5", "S-M010-1", "S-M010-2", "S-M010-3", "S-M010-4", "S-M010-5",
  "S-M010-6", "S-M009-1", "S-M009-2", "S-M011", "S-M012", "S-M005", "S-M006",
  "S-M004-4", "S-M004-5", "S-M004-6", "S-S013", "S-S014", "S-S015", "S-S016",
  "S-S017", "S-S018", "S-S019", "S-S020", "S-M007-1", "S-M007-2", "S-M007-3",
  "S-M007-4", "S-M007-5", "S-M007-6", "S-M008-1", "S-M008-2", "S-M009-3",
  "S-M009-4", "S-M009-5", "S-M009-6", "S-S021-1", "S-S021-2", "S-S021-3",
  "S-S021-4", "S-S021-5", "S-M018", "S-S029", "S-S030", "S-M023", "S-S031",
  "S-M024", "S-M017", "S-M019", "S-M020", "S-M021", "S-S028", "S-M022",
  "S-M001", "S-M030", "S-M031-1", "S-M031-2", "S-M032", "S-M033", "S-M034",
  "S-M035-1", "S-M035-2", "S-M036-1", "S-M036-2", "S-M037", "S-M038",
  "S-M039-1", "S-M039-2", "S-M040-1", "S-M040-2", "S-S001-2", "S-M069",
  "S-M070", "S-M071", "S-M072", "S-M073", "S-S077", "S-S078", "S-S079",
  "S-S080", "S-S081", "S-S082", "S-S083", "S-S084", "S-S037", "S-S038",
  "S-S039", "S-S040", "S-S041", "S-S042", "S-S043", "S-S044", "S-S045",
  "S-S046", "S-S047", "S-S048", "S-S049", "S-S050", "S-S051", "S-S052", "S-S053"
)
file_set <- patients_500k
merged_file <- "patients_500k.rds"
file_set <- patients_full
merged_file <- "patients_full.rds"
patient_metadata <- read.csv(meta_data_csv_file, header = TRUE,
                             row.names = 1)
cell_metadata <- read.csv("./cell_annotation.csv", header = TRUE,
                          row.names = 1)
# load BP cell matrices
data_list <- c()
metadata_list <- c()
for (i in (1:length(file_set))) {
  pat_id <- file_set[i]
  bp_dir <- paste0(files_dir, pat_id, "_data_BP")
  mat <- open_matrix_dir(dir = bp_dir)
  cmeta <- cell_metadata[colnames(mat), ]
  pmeta <- patient_metadata[pat_id, ]

  cmeta$pid <- pat_id
  cmeta$type <- pmeta[1, "characteristics..Sample.type"]
  cmeta$covid <- pmeta[1, "characteristics..SARS.CoV.2"]
  cmeta$batch_id <- pmeta[1, "characteristics...Datasets"]
  cmeta$severity <- pmeta[1, "characteristics..CoVID.19.severity"]
  cmeta$cell_type <- cmeta$celltype
  cmeta$cell_sub_type <- cmeta$celltype
  cmeta$cell_type <- cmeta$majorType

  metadata_list[[i]] <- cmeta
  data_list[[i]] <- mat
}
names(data_list) <- file_set
metadata <- Reduce(rbind, metadata_list)

# load BP cell matrices
options(Seurat.object.assay.version = "v5")
merged_object <- CreateSeuratObject(counts = data_list, meta.data = metadata)
merged_object

saveRDS(
  object = merged_object,
  file = merged_file,
  destdir = merged_dir
)
