import scanpy as sc
import pandas as pd
from pathlib import PurePath

# Generates the h5ad files from the input matrices
#
INPUT_DATA_DIR = "./aortic_valve"
OUTPUT_DATA_DIR = "./aortic_valve/"


def prep_dataset(data_dir, matrix_dir, meta_file,
                 batch_name, out_dir):
    meta_file = PurePath(data_dir, meta_file)
    matrix_dir = data_dir + "/" + matrix_dir
    meta_data = pd.read_csv(meta_file, sep="\t")  # type:ignore
    adx = sc.read_10x_mtx(matrix_dir)
    adx = adx[meta_data.cell_id, ]  # type:ignore
    meta_data['batch'] = batch_name  # type:ignore
    meta_data['cell_id'] = [batch_name + '-' + x for x in meta_data['cell_id']]  # type:ignore
    meta_data = meta_data.set_index(['cell_id'])  # type:ignore
    adx.obs = meta_data  # type:ignore
    adx.write(PurePath(out_dir, batch_name + ".h5ad"))


def main():
    prep_dataset(INPUT_DATA_DIR, "Healthy_1_filtered_feature_bc_matrix/",
                 "Healthy_1.tsv", "H1", OUTPUT_DATA_DIR)
    prep_dataset(INPUT_DATA_DIR, "Healthy_2_filtered_feature_bc_matrix/",
                 "Healthy_2.tsv", "H2", OUTPUT_DATA_DIR)
    prep_dataset(INPUT_DATA_DIR, "Calcified_3_filtered_feature_bc_matrix/",
                 "Calcified_3.tsv", "C3", OUTPUT_DATA_DIR)
    prep_dataset(INPUT_DATA_DIR, "Calcified_4_filtered_feature_bc_matrix/",
                 "Calcified_4.tsv", "C4", OUTPUT_DATA_DIR)


if __name__ == "__main__":
    main()
