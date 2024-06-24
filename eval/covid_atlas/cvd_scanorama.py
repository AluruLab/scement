import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
import scanorama as scanr
import argparse
import sys
import time
import pathlib
import scipy
from pathlib import PurePath


NGENES_MAP = {"500": 500, "1K": 1000, "5K": 5000, "10K": 10000, "full": None}
PATIENT_LIST = {
    "10K" : ['S-S062', 'S-M010-3', 'S-S054'],
    "25K" : ['S-S035-1', 'S-S037', 'S-S013', 'S-M077', 'S-S036-3'],
    "50K" : ['S-M018', 'S-S050', 'S-M066', 'S-S013', 'S-S035-1',
             'S-S001-2', 'S-S081', 'S-M037'],
    "150K": ['S-M061-2', 'S-M041-1', 'S-S021-4', 'S-M044-1', 'S-M043-2',
             'S-M053', 'S-S018', 'S-S036-3', 'S-S052', 'S-S035-3',
             'S-M056', 'S-S087-2', 'S-M049', 'S-M020', 'S-M001',
             'S-S082', 'S-M035-1', 'S-M012', 'S-S083', 'S-S050',
             'S-S027', 'S-M018', 'S-S086-2', 'S-S061'],
    "350K": [
        "S-S041", "S-S049", "S-M060-1", "S-S076-1", "S-M011", "S-M010-4", "S-S080",
        "S-M051", "S-S020", "S-S013", "S-S022-2", "S-S039", "S-M018", "S-M007-2",
        "S-M027", "S-M004-6", "S-M033", "S-M014", "S-S018", "S-S026", "S-S086-2",
        "S-S031", "S-M042-2", "S-S073-1", "S-M008-2", "S-S083", "S-S021-4",
        "S-S043", "S-M010-3", "S-S077", "S-M004-3", "S-M017", "S-S021-2",
        "S-M005", "S-M004-2", "S-M058-1", "S-S036-1", "S-M056", "S-S091",
        "S-S070-2", "S-M007-4", "S-M010-2", "S-M076-2", "S-M043-1",
        "S-M028", "S-S030", "S-S001-2", "S-S023", "S-S035-4", "S-M041-2",
        "S-M007-6", "S-S021-3", "S-S085-2", "S-S046", "S-M008-1",
        "S-S033", "S-M040-2", "S-M077", "S-S056", "S-M009-6"
    ],
    "500K": [
        "S-M043-1", "S-M071", "S-M004-3", "S-M079", "S-M007-1", "S-M055",
        "S-S087-2", "S-M014", "S-M077", "S-M051", "S-S015", "S-S035-4",
        "S-S047", "S-S085-2", "S-S049", "S-M026-1", "S-S040", "S-S083",
        "S-M072", "S-M053", "S-M035-2", "S-M026-3", "S-M009-4", "S-M016",
        "S-S034", "S-S063", "S-S051", "S-S031", "S-M048", "S-S025",
        "S-M004-6", "S-S026", "S-M043-2", "S-S054", "S-S036-3", "S-M015",
        "S-S084", "S-M066", "S-S090-2", "S-S020", "S-S080", "S-S022-3",
        "S-S021-4", "S-M059-1", "S-M007-3", "S-M007-5", "S-S077", "S-M004-2",
        "S-M009-5", "S-M017", "S-M063-1", "S-M044-1", "S-M037", "S-S023",
        "S-S070-3", "S-S076-2", "S-S013", "S-M010-4", "S-M059-2", "S-S024",
        "S-M041-1", "S-M039-1", "S-S027", "S-S043", "S-M040-1", "S-M026-2",
        "S-S037", "S-M067", "S-S028", "S-M004-5", "S-M042-2", "S-M036-1",
        "S-M062-2", "S-M004-1", "S-M021", "S-M022", "S-M042-1", "S-M039-2",
        "S-M018", "S-M078"
    ],
    "700K": [
        "S-S022-2", "S-S044", "S-S015", "S-M020", "S-M056", "S-M036-1",
        "S-S037", "S-S013", "S-S088-2", "S-M048", "S-M007-2", "S-M068",
        "S-M025", "S-S021-5", "S-S080", "S-S090-2", "S-M007-1", "S-S032-3",
        "S-S079", "S-S046", "S-S022-4", "S-S036-2", "S-M011", "S-M043-2",
        "S-M009-6", "S-M077", "S-M004-4", "S-S050", "S-M037", "S-S089-2",
        "S-M008-1", "S-S078", "S-M023", "S-M029", "S-S043", "S-M061-2",
        "S-S057", "S-M060-1", "S-M061-1", "S-M009-3", "S-S087-2", "S-M035-2",
        "S-S066", "S-S045", "S-M010-6", "S-S091", "S-S035-2", "S-M041-1",
        "S-S021-2", "S-M063-1", "S-M021", "S-S022-3", "S-S086-2", "S-M047",
        "S-M074-2", "S-S040", "S-S027", "S-M034", "S-M067", "S-S063",
        "S-S065", "S-S036-3", "S-M018", "S-S029", "S-M001", "S-S023",
        "S-S068", "S-M007-6", "S-M022", "S-S022-5", "S-S075-1", "S-M009-4",
        "S-M059-1", "S-S069-3", "S-S026", "S-S051", "S-M007-4", "S-M015",
        "S-S021-4", "S-M004-1", "S-M010-5", "S-M039-2", "S-M013", "S-S016",
        "S-M071", "S-S048", "S-S083", "S-M009-5", "S-M031-2", "S-M042-2",
        "S-M033", "S-S014", "S-M045", "S-M066", "S-M049", "S-M007-3",
        "S-S025", "S-M009-1", "S-M040-2", "S-S067", "S-S042", "S-M054",
        "S-M006", "S-M024", "S-M007-5", "S-S054", "S-S039", "S-M030",
        "S-M038", "S-M079", "S-S074-1", "S-S028", "S-S031", "S-M064", "S-S084"
    ],
    "full": [
        'S-S070-2', 'S-S070-3', 'S-S069-3', 'S-M056', 'S-M044-1', 'S-M043-1',
        'S-M048', 'S-M044-2', 'S-M043-2', 'S-S054', 'S-S056', 'S-M042-1',
        'S-M041-1', 'S-M049', 'S-M046', 'S-M047', 'S-S055', 'S-S057', 'S-M045',
        'S-M041-2', 'S-M042-2', 'S-M055', 'S-M053', 'S-S067', 'S-S065', 'S-M051',
        'S-S064', 'S-M054', 'S-M052', 'S-S068', 'S-S066', 'S-S059', 'S-S060',
        'S-M050', 'S-S061', 'S-S062', 'S-S063', 'S-M061-1', 'S-M061-2', 'S-S073-1',
        'S-S074-1', 'S-S074-2', 'S-M062-2', 'S-M058-1', 'S-M058-2', 'S-M063-1',
        'S-M063-2', 'S-M059-1', 'S-M059-2', 'S-M060-1', 'S-M060-2', 'S-S075-1',
        'S-S076-1', 'S-S076-2', 'S-S035-1', 'S-S035-2', 'S-S035-3', 'S-S035-4',
        'S-S036-1', 'S-S036-2', 'S-S036-3', 'S-M064', 'S-M066', 'S-M067', 'S-M068',
        'S-S091', 'S-S090-2', 'S-S092', 'S-M076-2', 'S-M074-2', 'S-S088-2',
        'S-S085-2', 'S-S089-2', 'S-S087-2', 'S-S086-2', 'S-M077', 'S-M078',
        'S-M079', 'S-S024', 'S-S026', 'S-M013', 'S-M014', 'S-M015', 'S-M016',
        'S-S023', 'S-S025', 'S-S027', 'S-S034', 'S-M025', 'S-S033', 'S-M028',
        'S-M029', 'S-M027', 'S-M026-1', 'S-S032-3', 'S-M026-2', 'S-M026-3',
        'S-M004-1', 'S-M004-2', 'S-M004-3', 'S-S022-1', 'S-S022-2', 'S-S022-3',
        'S-S022-4', 'S-S022-5', 'S-M010-1', 'S-M010-2', 'S-M010-3', 'S-M010-4',
        'S-M010-5', 'S-M010-6', 'S-M009-1', 'S-M009-2', 'S-M011', 'S-M012',
        'S-M005', 'S-M006', 'S-M004-4', 'S-M004-5', 'S-M004-6', 'S-S013', 'S-S014',
        'S-S015', 'S-S016', 'S-S017', 'S-S018', 'S-S019', 'S-S020', 'S-M007-1',
        'S-M007-2', 'S-M007-3', 'S-M007-4', 'S-M007-5', 'S-M007-6', 'S-M008-1',
        'S-M008-2', 'S-M009-3', 'S-M009-4', 'S-M009-5', 'S-M009-6', 'S-S021-1',
        'S-S021-2', 'S-S021-3', 'S-S021-4', 'S-S021-5', 'S-M018', 'S-S029',
        'S-S030', 'S-M023', 'S-S031', 'S-M024', 'S-M017', 'S-M019', 'S-M020',
        'S-M021', 'S-S028', 'S-M022', 'S-M001', 'S-M030', 'S-M031-1', 'S-M031-2',
        'S-M032', 'S-M033', 'S-M034', 'S-M035-1', 'S-M035-2', 'S-M036-1',
        'S-M036-2', 'S-M037', 'S-M038', 'S-M039-1', 'S-M039-2', 'S-M040-1',
        'S-M040-2', 'S-S001-2', 'S-M069', 'S-M070', 'S-M071', 'S-M072', 'S-M073',
        'S-S077', 'S-S078', 'S-S079', 'S-S080', 'S-S081', 'S-S082', 'S-S083',
        'S-S084', 'S-S037', 'S-S038', 'S-S039', 'S-S040', 'S-S041', 'S-S042',
        'S-S043', 'S-S044', 'S-S045', 'S-S046', 'S-S047', 'S-S048', 'S-S049',
        'S-S050', 'S-S051', 'S-S052', 'S-S053'
    ]
}

DATA_DIR = "/project/spc/i3/eval/covid_atlas/"
COVID_METADATA_XLSX = DATA_DIR + "/meta/covid_metadata.xlsx"
COVID_ANNOTATION_CSV = DATA_DIR + "/meta/cell_annotation.csv.gz"
GENE_MAPPING_10X_CSV = DATA_DIR + "/meta/gene_mapping.csv"


def load_covid_metadata(cvd_meta_file=COVID_METADATA_XLSX,
                        cvd_annot_file=COVID_ANNOTATION_CSV,
                        nmin=1000):
    pts_meta = pd.read_excel(cvd_meta_file, index_col="sampleID")  # type:ignore
    filter_types = ["fresh BALF", "fresh Sputum", "fresh PFMC"]
    pts_meta = pts_meta.loc[~pts_meta["Sample type"].isin(filter_types), ]
    # vcx = dfx["Sample type"].value_counts()
    #
    ptsqc_meta = pd.read_excel(cvd_meta_file, sheet_name="qc",  # type:ignore
                               index_col="sampleID")  # type:ignore
    ptsqc_meta = ptsqc_meta.loc[pts_meta.index]
    ptsqc_meta = ptsqc_meta.loc[ptsqc_meta["Ncells_manual_removal"] > nmin]
    pts_meta = pts_meta.loc[ptsqc_meta.index]
    cell_meta = pd.read_csv(cvd_annot_file, index_col='cellName')
    return ptsqc_meta, pts_meta, cell_meta


def filter_scale(adata, ntopg=5000):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if ntopg is not None:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=ntopg)


def load_covid_data_adlist(pts_meta, cell_meta, hd_dir, patients, ngenes):
    pts_ads = [ad.read_h5ad(hd_dir + "/" + x + ".h5ad") for x in patients]
    for i, x in enumerate(patients):
        pts_ads[i].obs["pid"] = x
        pts_ads[i].obs["batch"] = pts_meta.loc[x]["Batches"]
        pts_ads[i].obs["covid"] = pts_meta.loc[x]["SARS-CoV-2"]
        pts_ads[i].obs["severity"] = pts_meta.loc[x]["CoVID-19 severity"]
        pts_ads[i].obs["type"] = pts_meta.loc[x]['Sample type']
        pts_ads[i].obs["cell_type"] = cell_meta.loc[pts_ads[i].obs.index,
                                                    "majorType"]
        pts_ads[i].obs["cell_sub_type"] = cell_meta.loc[pts_ads[i].obs.index,
                                                        "celltype"]
        pts_ads[i].obs['source'] = 'COVID-STUDY'
        pts_ads[i].var.index = list(x for x in pts_ads[i].var["gene_id"])
        filter_scale(pts_ads[i], ngenes)
    return pts_ads


def load_covid_diseased_data_list(hd_dir, dis_kp, ngenes):
    _, dfx, cellfx = load_covid_metadata()
    if dis_kp is None:
        dis_kp = list(dfx[dfx['SARS-CoV-2'] == 'positive'].index)
    dis_ads = load_covid_data_adlist(dfx, cellfx, hd_dir, dis_kp, ngenes)
    return dis_ads


def scanorama_integrate2(scn_adatas, out_dir, file_pfx,
                         ntopg, plot_keys, workflow):
    tic = time.perf_counter()
    print(scn_adatas)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    scn_adatas_cor = scanr.correct_scanpy(scn_adatas, return_dimred=True)
    scanr.integrate_scanpy(scn_adatas_cor, dimred=50)
    toc = time.perf_counter()
    print(f"INTEGRATE in {toc - tic:0.4f} seconds")
    #
    tic = time.perf_counter()
    axcat = sc.concat(scn_adatas_cor, uns_merge="unique")
    axcat.obs_names_make_unique()
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    #
    if ntopg is not None:
        tic = time.perf_counter()
        sc.pp.highly_variable_genes(axcat, flavor="seurat",
                                    n_top_genes=ntopg)
        axcat = axcat[:, axcat.var.highly_variable]
        toc = time.perf_counter()
        print(f"HVG in {toc - tic:0.4f} seconds")
    #
    if workflow is True:
        tic = time.perf_counter()
        sc.pp.neighbors(axcat, use_rep="X_scanorama")
        sc.tl.umap(axcat)
        sc.tl.leiden(axcat, key_added="clusters")
        toc = time.perf_counter()
        print(f"POST PROC in {toc - tic:0.4f} seconds")
    else:
        plot_keys = None
    #
    tic = time.perf_counter()
    axcat_results_file = out_dir + "/scanorama_" + file_pfx + ".h5ad"
    axcat.write(PurePath(axcat_results_file))
    toc = time.perf_counter()
    print(f"SAVE-AXCAT in {toc - tic:0.4f} seconds")
    if plot_keys is None:
        return
    tic = time.perf_counter()
    axcat_plot_file = "scanorama_" + file_pfx + ".png"
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    sc.pl.umap(axcat, color=plot_keys, save=axcat_plot_file)
    toc = time.perf_counter()
    print(f"PLOT in {toc - tic:0.4f} seconds")
    return axcat


def combine_var_df(var_df_lst, col_name='gene_id',
                   update_cols=['n_cells', 'gene_name']):
    full_var_df = pd.merge(var_df_lst[0], var_df_lst[1], on=[col_name],
                           how='outer', suffixes=('_0', '_1'))
    enumlst = list(range(2, len(var_df_lst)))
    for ix in enumlst:
        full_var_df = pd.merge(full_var_df, var_df_lst[ix], on=[col_name],
                               how='outer', suffixes=('_0', '_1'))
        rename_map = {cx: cx + '_' + str(ix)
                      for cx in update_cols if cx in full_var_df.columns}
        if rename_map:
            full_var_df.rename(columns=rename_map,
                               inplace=True)
    full_var_df = full_var_df.fillna(0)
    return full_var_df


def gene_id2name_map(all_ncell_df, var_df_lst, col_name='gene_id'):
    gene_id2name = {x : [] for x in all_ncell_df[col_name]}
    for ix, vx in enumerate(var_df_lst):
        for _, y in vx[col_name].items():
            gene_id2name[y].append((y, ix))
        gids_set = set(x for x in vx[col_name])
        for y in gene_id2name.keys():
            if y not in gids_set:
                gene_id2name[y].append(('~', ix))
    return gene_id2name


def missing_genes_ad(adx, gene_id2name, col_name='gene_id'):
    adx_gset = set(adx.var[col_name])
    missing_ids = [x for x in gene_id2name.keys() if x not in adx_gset]
    missing_names = [gene_id2name[y][0][0] for y in missing_ids]
    missing_names = missing_ids
    obs_df = adx.obs.copy()
    var_df = pd.DataFrame({col_name: missing_ids})
    if 'n_cells' in adx.var.columns:
        var_df['n_cells'] = 0
    if 'gene_name' in adx.var.columns:
        var_df['gene_name'] = var_df[col_name]
    if 'highly_variable' in adx.var.columns:
        var_df['highly_variable'] = False
    var_df.index = missing_names
    # var_df.set_index(missing_names)
    fill_in_X = scipy.sparse.csr_matrix((obs_df.shape[0], len(missing_ids)),
                                        dtype=np.float32)
    fill_in_ad = ad.AnnData(fill_in_X, obs_df, var_df)
    return fill_in_ad


def fill_missing(adx, gene_id2name, col_name='gene_id'):
    aidx = [x for x in adx.var[col_name]]
    adx.var.index = aidx
    rdx = missing_genes_ad(adx, gene_id2name, col_name)
    conx = ad.concat([adx, rdx], axis=1)
    conx.obs = adx.obs
    return conx


def prep_combat_union(ad_objects):
    tic = time.perf_counter()
    au_objects = [ax[:, ~ax.var.gene_id.duplicated()] for ax in ad_objects]
    # Prepare the var data frame
    var_df_lst = [x.var[['gene_id', 'gene_name', 'n_cells', 'highly_variable']]
                  for x in au_objects]
    var_df_lst = [dfx.reset_index(drop=True) for dfx in var_df_lst]
    full_var_df = combine_var_df(var_df_lst, 'gene_id')
    gene_id2name = gene_id2name_map(full_var_df, var_df_lst, 'gene_id')
    srt_gene_id2name = {x: sorted(y) for x, y in gene_id2name.items()}
    gene_list = [x for x in gene_id2name.keys()]
    missing_ad_lst = [fill_missing(adx, srt_gene_id2name, "gene_id")
                      for adx in au_objects]
    missing_ad_lst2 = [x[:, gene_list] for x in missing_ad_lst]  # type:ignore
    # fill_var_df_lst = [x.var[['gene_id', 'n_cells', 'gene_name', 'highly_variable']]
    #                    for x in missing_ad_lst2]
    # final_var_df = combine_var_df(fill_var_df_lst, 'gene_id')
    # final_var_df.index = [x for x in final_var_df['gene_id']]  # type:ignore
    # axcat = ad.concat(missing_ad_lst2)
    # axcat.var = final_var_df
    # axcat.obs_names_make_unique()
    # print(au_objects)
    toc = time.perf_counter()
    print(f"UNION PREP in {toc - tic:0.4f} seconds")
    # print(axcat)
    return missing_ad_lst2


def main(in_args):
    data_dir = in_args.in_dir
    ngenes = NGENES_MAP[in_args.gxtype]
    dis_samples = PATIENT_LIST[in_args.cxtype]
    out_dir = in_args.out_dir
    out_prefix = "_".join(["run", in_args.cxtype, in_args.gxtype])
    #
    tic = time.perf_counter()
    au_objects = load_covid_diseased_data_list(data_dir, dis_samples, None)
    au_objects = [ax[:, ~ax.var.gene_id.duplicated()] for ax in au_objects]
    toc = time.perf_counter()
    print(f"LOAD in {toc - tic:0.4f} seconds")
    #
    if in_args.union:
        tic = time.perf_counter()
        var_df_lst = [x.var[['gene_id', 'n_cells']] for x in au_objects]
        full_var_df = combine_var_df(var_df_lst)
        gene_id2name = gene_id2name_map(full_var_df, var_df_lst)
        srt_gene_id2name = {x: sorted(y) for x, y in gene_id2name.items()}
        gene_list = [x for x in gene_id2name.keys()]
        missing_ad_lst = [fill_missing(adx, srt_gene_id2name) for adx in au_objects]
        missing_ad_lst2 = [x[:, gene_list] for x in missing_ad_lst]  # type:ignore
        # fill_var_df_lst = [x.var[['gene_id', 'n_cells']] for x in missing_ad_lst2]
        toc = time.perf_counter()
        print(f"FILL MISSING in {toc - tic:0.4f} seconds")
        au_objects = missing_ad_lst2
    print(au_objects)
    tic = time.perf_counter()
    plot_keys = None
    if in_args.plots:
        plot_keys = ['pid', 'covid', 'cell_type']

    axcat = scanorama_integrate2(au_objects, out_dir, out_prefix,
                                 ngenes, plot_keys, in_args.workflow)
    toc = time.perf_counter()
    print(f"COMPLETE WORKFLOW in {toc - tic:0.4f} seconds")
    print(axcat)


if __name__ == "__main__":
    rctypes = ["10K", "25K", "50K", "150K", "350K", "500K", "700K", "full"]
    rgtypes = list(NGENES_MAP.keys())
    default_in_dir = "./h5ad_full/"
    default_out_dir = "./scanorama_out/"
    parser = argparse.ArgumentParser(description='Integrate w. combat.')
    parser.add_argument("cxtype", choices=rctypes)
    parser.add_argument("gxtype", choices=rgtypes)
    parser.add_argument("-i", "--in_dir", default=default_in_dir)
    parser.add_argument("-o", "--out_dir", default=default_out_dir)
    parser.add_argument("-w", "--workflow", action="store_true")
    parser.add_argument("-l", "--plots", action="store_true")
    parser.add_argument("-u", "--union", action="store_true")
    #
    in_args = parser.parse_args(sys.argv[1:])
    print(in_args)
    main(in_args)
