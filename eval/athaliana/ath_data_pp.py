import scanpy as sc
import pandas as pd
import anndata as an
# from datetime import datetime
import sys
from pathlib import PurePath
#
# Datasets:
# 1. E-GEOD-152766
#     https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-152766/experiment-design
#    Shahan R, Hsu C, Nolan TM, Cole BJ, Taylor IW et al. (2020)
#     A single cell Arabidopsisroot atlas reveals developmental trajectories in wild type and cell identity mutants
# 2. E-GEOD-121619
#    https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-121619/downloads
#    Jean-Baptiste K, McFaline-Figueroa JL, Alexandre CM, Dorrity MW, Saunders L et al. (2019)
#     Dynamics of Gene Expression in Single Root Cells of Arabidopsis thaliana.
# 3. E-GEOD-123013
#    https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-123013/downloads
#    Ryu KH, Huang L, Kang HM, Schiefelbein J. (2019)
#      Single-Cell RNA Sequencing Resolves Molecular Relationships Among Individual Plant Cells.
# 4. E-GEOD-158761
#    https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-158761/downloads
#    Gala HP, Lanctot A, Jean-Baptiste K, Guiziou S, Chu JC et al. (2020)
#     A single cell view of the transcriptome during lateral root initiation inArabidopsis thaliana

meta_columns = ['Assay', 'access', 'ecotype', 'genotype', 'batch',
                'stress', 'age', 'stage', 'cell_type']


def get_common_genes(base_dir="./"):
    t66ad = an.read_h5ad(base_dir + "/E-GEOD-152766/E-GEOD-152766.h5ad")
    t19ad = an.read_h5ad(base_dir + "/E-GEOD-121619/E-GEOD-121619.h5ad")
    t13ad = an.read_h5ad(base_dir + "/E-GEOD-123013/E-GEOD-123013.h5ad")
    t61ad = an.read_h5ad(base_dir + "/E-GEOD-158761/E-GEOD-158761.h5ad")
    common_genes = list(set(t19ad.var_names) & set(t66ad.var_names) & set(t13ad.var_names) & set(t61ad.var_names))
    return common_genes


def filter_152766(common_genes, base_dir="./"):
    t66ad = an.read_h5ad(base_dir + "E-GEOD-152766/E-GEOD-152766.h5ad")
    t66_metax = t66ad.obs
    # t66 filters
    t66_label_filter = pd.isna(t66_metax['cell_type']) is False
    # cells of shr-2 genotype - 6942
    t66_shrgt_filter = ((t66_metax['genotype'] == 'shr-2') &
                        (t66_metax['batch'] == 'sc_52'))
    # cells of col-0 batch - 6396 cells
    t66_col0b_filter = (t66_metax['batch'] == 'col0') & t66_label_filter
    # cells of tnw1 batch - 4864 cells
    t66_tnw1b_filter = (t66_metax['batch'] == 'tnw1') & t66_label_filter
    # cells of sc_1 batch - 7097 cells
    t66_sc1b_filter = (t66_metax['batch'] == 'sc_1') & t66_label_filter
    #
    # t66 Filter 1
    t66_wt_filter = (t66_col0b_filter | t66_tnw1b_filter | t66_sc1b_filter)
    t66wtad = t66ad[t66_wt_filter, common_genes]
    t66wtad.write(PurePath(base_dir, "E-GEOD-152766/E-GEOD-152766-WT.h5ad"))
    # t66 Filter 2
    t66_batch_filter = (t66_col0b_filter | t66_shrgt_filter)
    t66btad = t66ad[t66_batch_filter, common_genes]
    t66btad.write(PurePath(base_dir, "E-GEOD-152766/E-GEOD-152766-BT.h5ad"))
    # t66 Filter 3 - atlas
    t66lblad = t66ad[t66_label_filter, common_genes]
    t66lblad.write(PurePath(base_dir, "E-GEOD-152766/E-GEOD-152766-LB.h5ad"))


def filter2_152766(base_dir="./"):
    t66ad = an.read_h5ad(base_dir + "E-GEOD-152766/E-GEOD-152766.h5ad")
    t66_metax = t66ad.obs
    # t66 filters
    t66_label_filter = pd.isna(t66_metax['cell_type']) is False
    t66_col0_filter = (t66_metax['ecotype'] == 'Col-0') & t66_label_filter
    print("Filter 2 : ", sum(t66_col0_filter))
    #
    # t66 Filter 3 - atlas
    t66col0 = t66ad[t66_col0_filter, ]  # type:ignore
    print(t66col0)
    t66col0.write(PurePath(base_dir, "E-GEOD-152766/E-GEOD-152766-COL0.h5ad"))


def prep_152766_gse(base_dir="./"):
    mtx_dir = base_dir + "GSE152766/"
    meta_file = PurePath(base_dir, "GSE152766/copilot_metadata.csv")
    t66_meta: pd.DataFrame = pd.read_csv(meta_file, dtype=object)  # type:ignore
    t66_meta['Assay'] = t66_meta['Cell_Barcode']
    t66_meta.index = t66_meta['Assay']  # type:ignore
    t66_meta['access'] = 'G152766'
    batches = list(set(t66_meta['orig.ident']) - set(['dc1', 'dc2', 'pp1']))  # type:ignore
    batch_meta_map = {x: t66_meta[t66_meta['orig.ident'] == x].copy() for x in batches}  # type:ignore
    batch_sfx_map = {x: list(set(z.split('_')[1] for z in dfx.Cell_Barcode))[0] for x, dfx in batch_meta_map.items()}
    batch_gt_map = {x: 'wild type genotype' for x in batches}
    batch_et_map = {x: 'Col-0' for x in batches}
    batch_gt_map['sc_25'] = 'scr-4'
    batch_gt_map['sc_36'] = 'scr-4'
    batch_gt_map['sc_52'] = 'shr-2'
    batch_gt_map['sc_53'] = 'shr-2'
    batch_et_map['sc_25'] = 'Landsberg x Col-0'
    batch_et_map['sc_36'] = 'Landsberg x Col-0'
    for bx in batches:
        batch_meta = batch_meta_map[bx]
        # batch_sfx = batch_sfx_map[bx]
        batch_meta['batch'] = bx
        batch_meta['stress'] = 'none'
        batch_meta['cell_type'] = batch_meta['celltype.anno']
        batch_meta['genotype'] = batch_gt_map[bx]
        batch_meta['ecotype'] = batch_et_map[bx]
        batch_meta['age'] = '5 day'
        batch_meta['stage'] = 'seedling developmental stage'
        batch_meta_map[bx] = batch_meta.loc[:, meta_columns].copy()

    batch_data_map = {x: sc.read_10x_mtx(mtx_dir + "/" + x + "/") for x in batches}
    for bx in batches:
        batch_data_map[bx].obs_names = [x + '_' + batch_sfx_map[bx] for x in batch_data_map[bx].obs_names]
    for bx in batches:
        batch_data_map[bx].obs = batch_meta_map[bx]  # type:ignore
    for bx in batches:
        out_file = PurePath(mtx_dir, bx + ".h5ad")
        batch_data_map[bx].write(out_file)  # "E-GEOD-152766/E-GEOD-152766.h5ad"


def prep_152766(base_dir="./"):
    #
    mtx_dir = base_dir + "E-GEOD-152766/data"
    meta_file = base_dir + "E-GEOD-152766/ExpDesign-E-GEOD-152766.tsv"
    out_file = base_dir + "E-GEOD-152766/E-GEOD-152766.h5ad"
    t66_data = sc.read_10x_mtx(mtx_dir)  # "E-GEOD-152766/data"
    t66_meta: pd.DataFrame = pd.read_csv(PurePath(meta_file),  # type:ignore
                                         sep="\t", dtype=object)  # type:ignore
    t66_meta.index = t66_meta['Assay']  # type:ignore
    t66_meta['access'] = 'G152766'
    t66_meta['genotype'] = t66_meta['Sample Characteristic[genotype]']
    t66_meta['ecotype'] = t66_meta['Sample Characteristic[ecotype]']
    t66_meta['stress'] = 'none'
    t66_meta['batch'] = t66_meta['access'] + '_' + t66_meta['Sample Characteristic[batch]']  # type:ignore
    # t66_gt_filter =  t66_meta['geno_type'] == 'wild type genotype'
    # t66_meta.loc[t66_gt_filter,'geno_type'] = 'col-0'
    t66_meta['age'] = t66_meta['Sample Characteristic[age]']
    t66_meta['stage'] = t66_meta['Sample Characteristic[developmental stage]']
    t66_meta['cell_type'] = t66_meta['Factor Value[inferred cell type - authors labels]']
    t66_metax = t66_meta.loc[:, meta_columns].copy()
    #
    t66_data.obs = t66_metax
    t66_data.write(PurePath(out_file))  # "E-GEOD-152766/E-GEOD-152766.h5ad"


def filter_121619(common_genes, base_dir="./"):
    t19ad = an.read_h5ad(base_dir + "E-GEOD-121619/E-GEOD-121619.h5ad")
    #
    t19_metax = t19ad.obs
    t19_label_filter = pd.isna(t19_metax['cell_type']) is False
    t19lblad = t19ad[t19_label_filter, common_genes]
    t19lblad.write(PurePath(base_dir, "E-GEOD-121619/E-GEOD-121619-LB.h5ad"))
    t19lbladfx = t19ad[t19_label_filter, common_genes]
    t19lbladfx.write(PurePath(base_dir, "E-GEOD-121619/E-GEOD-121619-FULL-LB.h5ad"))
    t19_cell_filter1 = (t19lblad.obs["cell_type"] != "root hair") & (t19lblad.obs["cell_type"] != "companion cell")
    t19lbladf1 = t19lblad[t19_cell_filter1, common_genes]
    t19lbladf1.write(PurePath(base_dir, "E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT1.h5ad"))


def filter2_121619(base_dir="./"):
    #
    t19ad = an.read_h5ad(base_dir + "E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT1.h5ad")
    t19ad = t19ad[t19ad.obs.batch != "G121619_heat", ]  # type:ignore
    t19ad.obs["batch"] = t19ad.obs["batch"].cat.rename_categories(
        {"G121619_none": "Batch-8"})  # type:ignore
    t19ad.write(PurePath(base_dir, "E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT2.h5ad"))


def prep_121619(base_dir="./"):
    mtx_dir = base_dir + "E-GEOD-121619/data"
    meta_file = base_dir + "E-GEOD-121619/ExpDesign-E-GEOD-121619.tsv"
    out_file = base_dir + "E-GEOD-121619/E-GEOD-121619.h5ad"
    t19_data = sc.read_10x_mtx(mtx_dir)  # "E-GEOD-121619/data")
    t19_meta: pd.DataFrame = pd.read_csv(meta_file,  # type:ignore
                                         sep="\t", dtype=object)  # type:ignore
    t19_meta.index = t19_meta['Assay']  # type:ignore
    t19_meta['access'] = 'G121619'
    t19_meta['genotype'] = t19_meta['Sample Characteristic[genotype]']
    t19_meta['ecotype'] = t19_meta['Sample Characteristic[ecotype]']
    t19_meta['stress'] = t19_meta['Factor Value[environmental stress]']
    t19_meta.loc[t19_meta['stress'] == 'warm/hot temperature regimen',
                 'stress'] = 'heat'
    t19_meta['batch'] = t19_meta['access'] + '_' + t19_meta['stress']  # type:ignore
    t19_meta['age'] = t19_meta['Sample Characteristic[age]']
    t19_meta['stage'] = t19_meta['Sample Characteristic[developmental stage]']
    t19_meta['cell_type'] = t19_meta['Factor Value[inferred cell type - ontology labels]']
    t19_metax = t19_meta.loc[:, meta_columns].copy()
    print(len(set(t19_metax['cell_type'])), set(t19_metax['cell_type']))
    t19_metax.loc[t19_metax['cell_type'] == 'columella root cap cell',
                  'cell_type'] = 'columella root cap'
    t19_metax.loc[t19_metax['cell_type'] == 'cortex cell', 'cell_type'] = 'cortex'
    t19_metax.loc[t19_metax['cell_type'] == 'endodermal cell',
                  'cell_type'] = 'endodermis'
    # t19_metax.loc[t19_metax['cell_type'] == 'non-hair root epidermal cell',
    #               'cell_type'] = 'non-hair root epidermis'
    t19_metax.loc[t19_metax['cell_type'] == 'non-hair root epidermal cell',
                  'cell_type'] = 'epidermis'
    t19_metax.loc[t19_metax['cell_type'] == 'phloem cell',
                  'cell_type'] = 'phloem'
    t19_metax.loc[t19_metax['cell_type'] == 'root hair cell',
                  'cell_type'] = 'root hair'
    t19_metax.loc[t19_metax['cell_type'] == 'xylem cell', 'cell_type'] = 'xylem'
    t19_metax.loc[t19_metax['cell_type'] == 'xylem pole pericycle cell',
                  'cell_type'] = 'xylem pole pericycle'
    t19_metax.loc[t19_metax['cell_type'] == 'phloem pole pericycle cell',
                  'cell_type'] = 'phloem pole pericycle'
    print(len(set(t19_metax['cell_type'])), set(t19_metax['cell_type']))
    t19_data.obs = t19_metax
    t19_data.write(PurePath(out_file))


def filter_123013(common_genes, base_dir="./"):
    t13ad = an.read_h5ad(base_dir + "E-GEOD-123013/E-GEOD-123013.h5ad")
    #
    t13_metax = t13ad.obs
    t13_label_filter = pd.isna(t13_metax['cell_type']) is False
    t13lblad = t13ad[t13_label_filter, common_genes]
    t13lblad.write(base_dir + "E-GEOD-123013/E-GEOD-123013-LB.h5ad")  # type:ignore
    #
    #
    t13_rhd6_filter = t13_metax['genotype'] == 'rhd6 mutant'
    t13rhdad = t13ad[t13_rhd6_filter, common_genes]
    t13rhdad.write(base_dir + "E-GEOD-123013/E-GEOD-123013-RHD.h5ad")  # type:ignore


def prep_123013(base_dir="./"):
    mtx_dir = base_dir + "E-GEOD-123013/data"
    meta_file = PurePath(base_dir, "E-GEOD-123013/ExpDesign-E-GEOD-123013.tsv")
    out_file = base_dir + "E-GEOD-123013/E-GEOD-123013.h5ad"
    t13_data = sc.read_10x_mtx(mtx_dir)
    t13_meta: pd.DataFrame = pd.read_csv(meta_file,
                                         sep="\t", dtype=object)  # type:ignore
    t13_meta.index = t13_meta['Assay']  # type:ignore
    t13_meta['access'] = 'G123013'
    t13_meta['genotype'] = t13_meta['Sample Characteristic[genotype]']
    t13_meta['ecotype'] = t13_meta['Sample Characteristic[ecotype]']
    t13_meta['stress'] = 'none'
    t13_meta['batch'] = t13_meta['access'] + '_' + t13_meta['genotype']  # type:ignore
    # t66_gt_filter =  t66_meta['geno_type'] == 'wild type genotype'
    # t66_meta.loc[t66_gt_filter,'geno_type'] = 'col-0'
    t13_meta['age'] = t13_meta['Sample Characteristic[age]']
    t13_meta['stage'] = t13_meta['Sample Characteristic[developmental stage]']
    t13_meta['cell_type'] = t13_meta['Factor Value[inferred cell type - authors labels]']
    t13_metax = t13_meta.loc[:, meta_columns].copy()
    print(len(set(t13_metax['cell_type'])), set(t13_metax['cell_type']))
    t13_metax.loc[t13_metax['cell_type'] == 'columella root cap cell 1',
                  'cell_type'] = 'columella root cap'
    t13_metax.loc[t13_metax['cell_type'] == 'non hair root epidermal cell 3',
                  'cell_type'] = 'non-hair root epidermis'
    t13_metax.loc[t13_metax['cell_type'] == 'pericycle 0',
                  'cell_type'] = 'pericycle'
    t13_metax.loc[t13_metax['cell_type'] == 'pericycle 5',
                  'cell_type'] = 'pericycle'
    t13_metax.loc[t13_metax['cell_type'] == 'root endodermis 2',
                  'cell_type'] = 'endodermis'
    t13_metax.loc[t13_metax['cell_type'] == 'root endodermis 7',
                  'cell_type'] = 'endodermis'
    t13_metax.loc[t13_metax['cell_type'] == 'root epidermal cell 8',
                  'cell_type'] = 'epidermis'
    t13_metax.loc[t13_metax['cell_type'] == 'root hair cell 4',
                  'cell_type'] = 'root hair'
    t13_metax.loc[t13_metax['cell_type'] == 'root cortex endodermis initial cell 6',
                  'cell_type'] = 'cortex endodermis initial'
    print(len(set(t13_metax['cell_type'])), set(t13_metax['cell_type']))
    t13_data.obs = t13_metax
    t13_data.write(PurePath(out_file))  # "E-GEOD-123013/E-GEOD-123013.h5ad")


def filter_158761(common_genes, base_dir="./"):
    #
    t61ad = an.read_h5ad(base_dir + "E-GEOD-158761/E-GEOD-158761.h5ad")
    t61_metax = t61ad.obs
    t61_label_filter = pd.isna(t61_metax['cell_type']) is False
    t61lblad = t61ad[t61_label_filter, common_genes]
    t61lblad.write(PurePath(base_dir, "E-GEOD-158761/E-GEOD-158761-LB.h5ad"))


def prep_158761_gse(base_dir="./"):
    mtx_dir = base_dir + "/GSE158761/mtx"
    meta_file = PurePath(base_dir, "GSE158761/barcodes_annot.tsv")
    out_file = PurePath(base_dir, "GSE158761/GSE158761.h5ad")
    fout_file1 = PurePath(base_dir, "GSE158761/GSE158761-FLT1.h5ad")
    t61_data: an.AnnData = sc.read_10x_mtx(mtx_dir)
    t61_meta: pd.DataFrame = pd.read_csv(meta_file,
                                         sep="\t", dtype=object)  # type:ignore
    t61_meta['Assay'] = t61_meta['barcode']
    t61_meta.index = t61_meta['Assay']  # type:ignore
    t61_meta['access'] = 'G158761'
    t61_meta['batch'] = ['Batch-' + x.split('-')[1] for x in t61_meta['barcode']]  # type:ignore
    t61_meta['stress'] = 'none'
    t61_meta['genotype'] = 'wild type genotype'
    t61_meta['ecotype'] = 'Col-0'
    t61_meta['age'] = '2 day'
    t61_meta['stage'] = 'seedling developmental stage'
    t61_meta['age'] = t61_meta['age'] + '_' + t61_meta['time']  # type:ignore
    t61_metax = t61_meta.loc[:, meta_columns].copy()
    t61_metax.loc[t61_metax['cell_type'] == 'Epidermis',
                  'cell_type'] = 'epidermis'
    t61_metax.loc[t61_metax['cell_type'] == 'Phloem Pole Pericycle',
                  'cell_type'] = 'phloem pole pericycle'
    t61_metax.loc[t61_metax['cell_type'] == 'Columella/ Root Cap',
                  'cell_type'] = 'columella root cap'
    t61_metax.loc[t61_metax['cell_type'] == 'Cortex',
                  'cell_type'] = 'cortex'
    t61_metax.loc[t61_metax['cell_type'] == 'Endodermis',
                  'cell_type'] = 'endodermis'
    t61_metax.loc[t61_metax['cell_type'] == 'Phloem',
                  'cell_type'] = 'phloem'
    t61_metax.loc[t61_metax['cell_type'] == 'Xylem',
                  'cell_type'] = 'xylem'
    t61_metax.loc[t61_metax['cell_type'] == 'Mature Pericycle',
                  'cell_type'] = 'pericycle'
    t61_metax.loc[t61_metax['cell_type'] == 'Xylem Pole Pericycle',
                  'cell_type'] = 'xylem pole pericycle'
    t61_metax.loc[t61_metax['cell_type'] == 'Lateral Root Endodermis',
                  'cell_type'] = 'endodermis'
    t61_metax.loc[t61_metax['cell_type'] == 'Lateral Root Primordia',
                  'cell_type'] = 'lateral root primordia'
    t61_metax.loc[t61_metax['cell_type'] == 'Ambiguous Stele Cells',
                  'cell_type'] = 'ambig. stele cells'

    t61_data.obs = t61_metax
    t61_data.write(out_file)

    filtered_t61_data = t61_data[((t61_data.obs["cell_type"] != "ambig. stele cells") &
                                  (t61_data.obs["cell_type"] != "lateral root primordia"))]
    filtered_t61_data.write(fout_file1)


def prep_158761(base_dir="./"):
    mtx_dir = base_dir + "E-GEOD-158761/data"
    meta_file = PurePath(base_dir, "E-GEOD-158761/ExpDesign-E-GEOD-158761.tsv")
    out_file = PurePath(base_dir, "E-GEOD-158761/E-GEOD-158761.h5ad")

    t61_data: an.AnnData = sc.read_10x_mtx(mtx_dir)
    t61_meta: pd.DataFrame = pd.read_csv(meta_file, sep="\t",
                                         dtype=object)  # type:ignore
    t61_meta.index = t61_meta['Assay']  # type:ignore
    t61_meta['access'] = 'G158761'
    t61_meta['genotype'] = t61_meta['Sample Characteristic[genotype]']
    t61_meta['ecotype'] = t61_meta['Sample Characteristic[ecotype]']
    t61_meta['stress'] = 'none'
    t61_meta['batch'] = t61_meta['access'] + '_' + t61_meta['genotype']  # type:ignore
    # t66_gt_filter =  t66_meta['geno_type'] == 'wild type genotype'
    # t66_meta.loc[t66_gt_filter,'geno_type'] = 'col-0'
    t61_meta['age'] = t61_meta['Sample Characteristic[age]']
    t61_meta['stage'] = t61_meta['Sample Characteristic[developmental stage]']
    t61_meta['cell_type'] = t61_meta['Factor Value[inferred cell type - authors labels]']
    t61_metax = t61_meta.loc[:, meta_columns].copy()
    print(len(set(t61_metax['cell_type'])), set(t61_metax['cell_type']))
    t61_metax.loc[t61_metax['cell_type'] == 'columella/ root cap',
                  'cell_type'] = 'columella root cap'
    print(len(set(t61_metax['cell_type'])), set(t61_metax['cell_type']))
    t61_data.obs = t61_metax
    t61_data.write(out_file)


def main(rtype, in_dset, base_dir):
    if rtype == "prep":
        if in_dset == "152766":
            prep_152766(base_dir)
        elif in_dset == "121619":
            prep_121619(base_dir)
        elif in_dset == "123013":
            prep_123013(base_dir)
        elif in_dset == "158761":
            prep_158761(base_dir)
    elif rtype == "prep2":
        if in_dset == "152766":
            prep_152766_gse(base_dir)
        elif in_dset == "158761":
            prep_158761_gse(base_dir)
    elif rtype == "filter":
        common_genes = get_common_genes()
        if in_dset == "152766":
            filter_152766(common_genes, base_dir)
        elif in_dset == "152766":
            filter_121619(common_genes, base_dir)
        elif in_dset == "123013":
            filter_123013(common_genes, base_dir)
        elif in_dset == "158761":
            filter_158761(common_genes, base_dir)
    elif rtype == "filter2":
        if in_dset == "152766":
            filter2_152766(base_dir)
        if in_dset == "121619":
            filter2_121619(base_dir)


if __name__ == "__main__":
    run_types = ["prep", "filter", "filter2"]
    dsets = ["152766", "121619", "123013", "158761"]
    usage_str = "Usage : " + sys.argv[0] + " "
    usage_str += "/".join(run_types) + " "
    usage_str += "/".join(dsets) + " "
    usage_str += "base_dir"
    if len(sys.argv) <= 3:
        print(usage_str)
    else:
        rtype = sys.argv[1]
        in_dset = sys.argv[2]
        base_dir = sys.argv[3]
        if (in_dset not in dsets or rtype not in run_types):
            print(usage_str)
        else:
            main(rtype, in_dset, base_dir)
