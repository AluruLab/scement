import anndata as ad
import numpy as np
import pandas as pd
import argparse
import json
import scipy
import time


def combine_var_df(var_df_lst, col_name='gene_ids'):
    full_var_df = pd.merge(var_df_lst[0], var_df_lst[1], on=[col_name],
                           how='outer', suffixes=('_0', '_1'))
    enumlst = list(range(2, len(var_df_lst)))
    for ix in enumlst:
        full_var_df = pd.merge(full_var_df, var_df_lst[ix], on=[col_name],
                               how='outer', suffixes=('_0', '_1'))
        if 'n_cells' in full_var_df.columns:
            full_var_df.rename(columns={'n_cells': 'n_cells_' + str(ix)},
                               inplace=True)
    full_var_df = full_var_df.fillna(0)
    return full_var_df


def gene_id2name_map(all_ncell_df, var_df_lst, col_name='gene_ids'):
    gene_id2name = {x : [] for x in all_ncell_df[col_name]}
    for ix, vx in enumerate(var_df_lst):
        for x, y in vx[col_name].items():
            gene_id2name[y].append((x, ix))
        gids_set = set(x for x in vx[col_name])
        for y in gene_id2name.keys():
            if y not in gids_set:
                gene_id2name[y].append(('~', ix))
    return gene_id2name


def missing_genes_ad(adx, gene_id2name, col_name='gene_ids'):
    adx_gset = set(adx.var[col_name])
    adx_gset = set(adx.var[col_name])
    missing_ids = [x for x in gene_id2name.keys() if x not in adx_gset]
    missing_names = [gene_id2name[y][0][0] for y in missing_ids]
    missing_names = missing_ids
    obs_df = adx.obs.copy()
    var_df = pd.DataFrame({col_name: missing_ids})
    if 'n_cells' in var_df.columns:
        var_df['n_cells'] = 0
    var_df.index = missing_names
    fill_in_X = scipy.sparse.csr_matrix((obs_df.shape[0], len(missing_ids)),
                                        dtype=np.float32)
    fill_in_ad = ad.AnnData(fill_in_X, obs_df, var_df)
    return fill_in_ad


def fill_missing(adx, gene_id2name, col_name='gene_ids'):
    aidx = [x for x in adx.var[col_name]]
    adx.var.index = aidx
    rdx = missing_genes_ad(adx, gene_id2name)
    conx = ad.concat([adx, rdx], axis=1)
    conx.obs = adx.obs
    return conx


def main(json_data):
    data_dir = json_data['DATADIR']
    in_files = json_data['ADFILE']
    # batch_names = json_data['BATCH']
    out_dir = json_data['OUTDIR']
    out_file = json_data['OUTFILE']
    results_file = out_dir + "/" + out_file
    #
    #
    tic = time.perf_counter()
    ad_objects = [ad.read_h5ad(data_dir + "/" + fx) for fx in in_files]
    # var_lst = [x.var for x in ad_objects]
    var_df_lst = [x.var[['gene_ids', 'n_cells']] for x in ad_objects]
    full_var_df = combine_var_df(var_df_lst)
    toc = time.perf_counter()
    print(f"LOAD in {toc - tic:0.4f} seconds")
    #
    #
    tic = time.perf_counter()
    gene_id2name = gene_id2name_map(full_var_df, var_df_lst)
    srt_gene_id2name = {x: sorted(y) for x, y in gene_id2name.items()}
    gene_list = [x for x in gene_id2name.keys()]
    missing_ad_lst = [fill_missing(adx, srt_gene_id2name) for adx in ad_objects]
    missing_ad_lst2 = [x[:, gene_list] for x in missing_ad_lst]  # type:ignore
    fill_var_df_lst = [x.var[['gene_ids', 'n_cells']] for x in missing_ad_lst2]
    final_var_df = combine_var_df(fill_var_df_lst)
    final_var_df.index = [x for x in final_var_df['gene_ids']]
    final_conx2 = ad.concat(missing_ad_lst2)
    final_conx2.obs_names_make_unique()
    final_conx2.var = final_var_df
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(final_conx2)
    #
    #
    tic = time.perf_counter()
    final_conx2.write(results_file)
    toc = time.perf_counter()
    print(f"SAVE in {toc - tic:0.4f} seconds")


if __name__ == "__main__":
    help_str = """Path to input json file, File expected to have following format:
{
   "DATADIR": "/storage/home/hcoda1/0/spc/data2/",
  "OUTDIR" : "/storage/home/hcoda1/0/spc/data2/scr_data/rna_combo/",
   "ADFILE" : [
       "scr_data/rna_3p_chrx_10k_n02/pre_processed.h5ad",
       "scr_data/rna_3p_chrx_ht_20K_n03/pre_processed.h5ad",
       "scr_data/rna_3p_ctrlr_10k_n01/pre_processed.h5ad",
       "scr_data/rna_chrom_conn1_5k_n07/pre_processed.h5ad",
       "scr_data/rna_chrom_conn5_5k_n06/pre_processed.h5ad",
       "scr_data/rna_chromv2_10k_n08/pre_processed.h5ad",
       "scr_data/rna_chromv2_v2chem_1k_n09/pre_processed.h5ad",
       "scr_data/rna_chromv2_v3chem_1k_n10/pre_processed.h5ad",
       "scr_data/rna_cv1_68K_n11/pre_processed.h5ad",
       "scr_data/rna_cv1_donora_n12/pre_processed.h5ad",
       "scr_data/rna_cv1_donorb_n13/pre_processed.h5ad",
       "scr_data/rna_cv1_donorc_n14/pre_processed.h5ad",
       "scr_data/rna_tgt_10k_n05/pre_processed.h5ad",
       "scr_data/rna_pbmc_600k/pre_processed.h5ad"
   ],
  "BATCH" :  ["rna_3p_chrx_10k_n02", "rna_3p_chrx_ht_20K_n03",
       "rna_3p_ctrlr_10k_n01", "rna_chrom_conn1_5k_n07",
       "rna_chrom_conn5_5k_n06", "rna_chromv2_10k_n08",
       "rna_chromv2_v2chem_1k_n09", "rna_chromv2_v3chem_1k_n10",
       "rna_cv1_68K_n11", "rna_cv1_donora_n12", "rna_cv1_donorb_n13",
       "rna_cv1_donorc_n14", "rna_tgt_10k_n05", "rna_pbmc_600k" ],
  "OUTFILE" : "combo_joint_matrix.h5ad"
}
            """
    parser = argparse.ArgumentParser(description='Combine datasets with union.')
    parser.add_argument('json_file', type=str, help=help_str)
    args = parser.parse_args()
    with open(args.json_file) as f:
        json_data = json.load(f)
    main(json_data)
