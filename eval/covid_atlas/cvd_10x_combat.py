import pandas as pd
import scanpy as sc
import anndata as an
import numpy as np
import scipy
import scipy.sparse
import sys
import time
import pathlib
import argparse
import select_data as sd
import scement as sct

DATA_DIR = "/project/spc/i3/covid-scrnaseq/"
COVID_METADATA_XLSX = DATA_DIR + "/meta/covid_metadata.xlsx"
COVID_ANNOTATION_CSV = DATA_DIR + "/meta/cell_annotation.csv.gz"
GENE_MAPPING_10X_CSV = DATA_DIR + "/meta/gene_mapping.csv"
COMBO_10X_H5AD = DATA_DIR + "/data/combo_intersect.h5ad"
COMBO4_10X_H5AD = DATA_DIR + "/data/combo_intersect4.h5ad"
COMBO4_10X_META_CSV = DATA_DIR + "/data/combo_intersect4_mdata.csv"
COMBO7_10X_H5AD = DATA_DIR + "/data/combo_intersect7.h5ad"
COMBO7_10X_META_CSV = DATA_DIR + "/data/combo_intersect7_mdata.csv"


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


def load_10X_350K_data(full_ad_file=COMBO_10X_H5AD,
                       gm_file=GENE_MAPPING_10X_CSV):
    combo4_batches = [
        'rna_chromv2_10k_n08', 'rna_cv1_donorb_n13', 'rna_3p_ctrlr_10k_n01',
        'rna_chrom_conn1_5k_n07', 'rna_3p_chrx_ht_20K_n03',
        'rna_chrom_conn5_5k_n06', 'rna_tgt_10k_n05', 'rna_3p_chrx_10k_n02',
        'rna_cv1_donora_n12']
    tic = time.perf_counter()
    adx = an.read_h5ad(full_ad_file)
    lkdf = pd.read_csv(gm_file, index_col=0)
    new_obs = adx.obs
    new_obs['pid'] = new_obs['batch']
    adx.obs = new_obs
    obs_600k = new_obs.loc[new_obs["batch"] == "rna_pbmc_600k"]
    sample_600k = list(obs_600k.sample(n=265000, replace=False).index)
    left_xk = list(new_obs.loc[new_obs["batch"].isin(combo4_batches)].index)
    sel_lsx = left_xk + sample_600k
    adx = adx[sel_lsx]  # type:ignore
    lkdf2 = lkdf.copy()
    lkdf2.index = lkdf2.gene_id
    # lkdf2.set_index("gene_id", drop=False)
    new_var = lkdf2.loc[adx.var.index, ]
    new_var.index = new_var.gene_name
    new_var['gene_id2'] = new_var.gene_id
    new_var['gene_id'] = new_var.gene_name
    adx.var = new_var
    toc = time.perf_counter()
    print(f"10X 50K LOAD in {toc - tic:0.4f} seconds")
    return adx


def load_10X_50K_data(ad_file=COMBO4_10X_H5AD,
                      md_file=COMBO4_10X_META_CSV,
                      gm_file=GENE_MAPPING_10X_CSV):
    tic = time.perf_counter()
    adx = an.read_h5ad(ad_file)
    mdx = pd.read_csv(md_file)
    cxnames = list(mdx.columns)
    cxnames[0] = 'cell_id'
    cxmap = {x: y for x, y in zip(list(mdx.columns), cxnames)}
    mdx = mdx.rename(columns=cxmap)
    mdx.index = mdx.cell_id
    # mdx.set_index("cell_id", drop=False)
    mdx['pid'] = mdx.batch
    mdx['source'] = '10X-REPO'
    mdx = mdx.drop('batch', axis=1)
    new_obs = pd.concat([adx.obs, mdx], axis=1)
    new_obs['cell_type'] = new_obs['predicted.celltype.l1']
    adx.obs = new_obs
    lkdf = pd.read_csv(gm_file, index_col=0)
    lkdf2 = lkdf.copy()
    lkdf2.index = lkdf2.gene_id
    # lkdf2.set_index("gene_id", drop=False)
    new_var = lkdf2.loc[adx.var.index, ]
    new_var.index = new_var.gene_name
    new_var['gene_id2'] = new_var.gene_id
    new_var['gene_id'] = new_var.gene_name
    adx.var = new_var
    toc = time.perf_counter()
    print(f"10X 50K LOAD in {toc - tic:0.4f} seconds")
    return adx


def load_10X_180K_data(
        ad_file=COMBO7_10X_H5AD,
        md_file=COMBO7_10X_META_CSV,
        gm_file=GENE_MAPPING_10X_CSV):
    tic = time.perf_counter()
    adx = an.read_h5ad(ad_file)
    mdx = pd.read_csv(md_file)
    cxnames = list(mdx.columns)
    cxnames[0] = 'cell_id'
    cxmap = {x: y for x, y in zip(list(mdx.columns), cxnames)}
    mdx = mdx.rename(columns=cxmap)
    mdx.index = mdx.cell_id
    # mdx.set_index("cell_id", drop=False)
    mdx['pid'] = mdx.batch
    mdx['source'] = '10X-REPO'
    mdx = mdx.drop('batch', axis=1)
    new_obs = pd.concat([adx.obs, mdx], axis=1)
    new_obs['cell_type'] = new_obs['predicted.celltype.l1']
    adx.obs = new_obs
    lkdf = pd.read_csv(gm_file, index_col=0)
    lkdf2 = lkdf.copy()
    lkdf2.index = lkdf2.gene_id
    # lkdf2.set_index("gene_id", drop=False)
    new_var = lkdf2.loc[adx.var.index, ]
    new_var.index = new_var.gene_name
    new_var['gene_id2'] = new_var.gene_id
    new_var['gene_id'] = new_var.gene_name
    adx.var = new_var
    toc = time.perf_counter()
    print(f"10X 180K LOAD in {toc - tic:0.4f} seconds")
    return adx


def filter_scale(adata, ntopg=5000):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if ntopg is not None:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=ntopg)


def prep_combat_intersect(ad_objects):
    tic = time.perf_counter()
    au_objects = [ax[:, ~ax.var.gene_id.duplicated()] for ax in ad_objects]
    print(au_objects)
    axcat = an.concat(au_objects, merge='same')
    axcat.obs_names_make_unique()
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    return axcat


def combine_var_df(var_df_lst, col_name='gene_ids',
                   update_cols=['n_cells', 'gene_name', 'highly_variable']):
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


def gene_id2name_map(all_ncell_df, var_df_lst, col_name='gene_ids'):
    gene_id2name = {x : [] for x in all_ncell_df[col_name]}
    for ix, vx in enumerate(var_df_lst):
        for _, y in vx[col_name].items():
            gene_id2name[y].append((y, ix))
        gids_set = set(x for x in vx[col_name])
        for y in gene_id2name.keys():
            if y not in gids_set:
                gene_id2name[y].append(('~', ix))
    return gene_id2name


def missing_genes_ad(adx, gene_id2name, col_name='gene_ids'):
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
    fill_in_ad = an.AnnData(fill_in_X, obs_df, var_df)
    return fill_in_ad


def fill_missing(adx, gene_id2name, col_name='gene_ids'):
    aidx = [x for x in adx.var[col_name]]
    adx.var.index = aidx
    rdx = missing_genes_ad(adx, gene_id2name, col_name)
    conx = an.concat([adx, rdx], axis=1)
    conx.obs = adx.obs
    return conx


def prep_combat_union(ad_objects):
    tic = time.perf_counter()
    au_objects = [ax[:, ~ax.var.gene_id.duplicated()] for ax in ad_objects]
    # Prepare the var data frame
    var_sel_lst = ['gene_id', 'gene_name', 'n_cells', 'highly_variable']
    var_sel_lst = ['gene_id', 'gene_name', 'n_cells']
    var_df_lst = [x.var[var_sel_lst]
                  for x in au_objects]
    var_df_lst = [dfx.reset_index(drop=True) for dfx in var_df_lst]
    full_var_df = combine_var_df(var_df_lst, 'gene_id')
    gene_id2name = gene_id2name_map(full_var_df, var_df_lst, 'gene_id')
    srt_gene_id2name = {x: sorted(y) for x, y in gene_id2name.items()}
    gene_list = [x for x in gene_id2name.keys()]
    missing_ad_lst = [fill_missing(adx, srt_gene_id2name, "gene_id")
                      for adx in au_objects]
    missing_ad_lst2 = [x[:, gene_list] for x in missing_ad_lst]  # type:ignore
    fill_var_df_lst = [x.var[var_sel_lst]
                       for x in missing_ad_lst2]
    final_var_df = combine_var_df(fill_var_df_lst, 'gene_id')
    final_var_df.index = [x for x in final_var_df['gene_id']]  # type:ignore
    axcat = an.concat(missing_ad_lst2)
    axcat.var = final_var_df
    axcat.obs_names_make_unique()
    print(au_objects)
    toc = time.perf_counter()
    print(f"UNION in {toc - tic:0.4f} seconds")
    print(axcat)
    return axcat


def load_covid_data_adlist(pts_meta, cell_meta, hd_dir, patients, ngenes):
    pts_ads = [an.read_h5ad(hd_dir + "/" + x + ".h5ad") for x in patients]
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
        filter_scale(pts_ads[i], ngenes)
    return pts_ads


def load_covid_control_data(hd_dir, ngenes):
    tic = time.perf_counter()
    _, dfx, cellfx = load_covid_metadata()
    ctrl_kp = list(dfx[dfx['SARS-CoV-2'] == 'negative'].index)
    ctrl_ads = load_covid_data_adlist(dfx, cellfx, hd_dir, ctrl_kp, ngenes)
    axcat = prep_combat_intersect(ctrl_ads)
    toc = time.perf_counter()
    print(f"CTRL LOAD in {toc - tic:0.4f} seconds")
    return ctrl_ads, axcat


def load_covid_infected_data_list(hd_dir, dis_kp, ngenes):
    _, dfx, cellfx = load_covid_metadata()
    if dis_kp is None:
        dis_kp = list(dfx[dfx['SARS-CoV-2'] == 'positive'].index)
    dis_ads = load_covid_data_adlist(dfx, cellfx, hd_dir, dis_kp, ngenes)
    return dis_ads


def load_covid_full_data_list(hd_dir, dis_kp, ngenes):
    _, dfx, cellfx = load_covid_metadata()
    dis_kp = list(dfx.index)
    full_ads = load_covid_data_adlist(dfx, cellfx, hd_dir, dis_kp, ngenes)
    return full_ads


def load_covid_infected_data(hd_dir, dis_kp, ngenes):
    tic = time.perf_counter()
    dis_ads = load_covid_infected_data_list(hd_dir, dis_kp, ngenes)
    toc = time.perf_counter()
    print(f"DIS LOAD in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    axcat = prep_combat_intersect(dis_ads)
    toc = time.perf_counter()
    print(f"INTERSECT LOAD in {toc - tic:0.4f} seconds")
    return dis_ads, axcat


def load_covid_infected_data_union(hd_dir, dis_kp, ngenes):
    tic = time.perf_counter()
    dis_ads = load_covid_infected_data_list(hd_dir, dis_kp, ngenes)
    toc = time.perf_counter()
    print(f"DIS LOAD {ngenes} in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    axcat = prep_combat_union(dis_ads)
    toc = time.perf_counter()
    print(f"UNION LOAD in {toc - tic:1.4f} seconds")
    return dis_ads, axcat


def load_covid_full_data(hd_dir, dis_kp, ngenes):
    tic = time.perf_counter()
    full_ads = load_covid_full_data_list(hd_dir, dis_kp, ngenes)
    toc = time.perf_counter()
    print(f"DIS LOAD in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    axcat = prep_combat_intersect(full_ads)
    toc = time.perf_counter()
    print(f"INTERSECT LOAD in {toc - tic:0.4f} seconds")
    return full_ads, axcat


def load_covid_full_data_union(hd_dir, dis_kp, ngenes):
    tic = time.perf_counter()
    dis_ads = load_covid_full_data_list(hd_dir, dis_kp, ngenes)
    toc = time.perf_counter()
    print(f"DIS LOAD {ngenes} in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    axcat = prep_combat_union(dis_ads)
    toc = time.perf_counter()
    print(f"UNION LOAD in {toc - tic:1.4f} seconds")
    return dis_ads, axcat


def combat_and_plot(axcat, axcat_plot_file, axcat_file,
                    ntopg=5000, scement=True,
                    batch_key="pid", cvats=['batch', 'source'],
                    plot_keys=['pid', 'covid', 'cell_type'], post=False,
                    workflow=True, filter_pct=0.0):
    if workflow is False:
        plot_keys = None
    if (axcat is not None) and (filter_pct > 0.0):
        tic = time.perf_counter()
        pre_batches = set(axcat.obs[batch_key])
        celldx, genedx = sc.pp.calculate_qc_metrics(axcat)  # type:ignore
        print(celldx)
        print(genedx)
        npct = filter_pct
        ngenes = genedx.shape[0]
        nzcells_flag = celldx.n_genes_by_counts > int(ngenes * npct)
        nzpct = (100 * float(sum(nzcells_flag))) / celldx.shape[0]
        axcat = axcat[nzcells_flag, :]
        sc.pp.filter_genes(axcat, min_cells=3)
        post_batches = set(axcat.obs[batch_key])
        toc = time.perf_counter()
        print("NPCT B", pre_batches, post_batches, pre_batches - post_batches, nzpct)
        print("NPCT 2", axcat)
        print(f"PRE-FILTER in {toc - tic:0.4f} seconds")
    #
    if (post is False) and (ntopg is not None) and (ntopg > 0):
        tic = time.perf_counter()
        sc.pp.highly_variable_genes(axcat, flavor="seurat", n_top_genes=ntopg)
        axcat = axcat[:, axcat.var.highly_variable]
        toc = time.perf_counter()
        print(f"HVG in {toc - tic:0.4f} seconds")
    print(axcat)
    tic = time.perf_counter()
    if scement:
        sct.sct_sparse(axcat, key=batch_key, covariates=cvats, inplace=True)
    else:
        sct.combat(axcat, key=batch_key, covariates=cvats, inplace=True)
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")
    if (post is True) and (ntopg is not None) and (ntopg > 0):
        tic = time.perf_counter()
        sc.pp.highly_variable_genes(axcat, flavor="cell_ranger", n_top_genes=ntopg)
        axcat = axcat[:, axcat.var.highly_variable]
        toc = time.perf_counter()
        print(f"HVG in {toc - tic:0.4f} seconds")
    if workflow:
        # Main Workflow
        tic = time.perf_counter()
        sc.tl.pca(axcat, svd_solver='arpack', n_comps=20)
        toc = time.perf_counter()
        print(f"PCA in {toc - tic:0.4f} seconds")
        #
        tic = time.perf_counter()
        sc.pp.neighbors(axcat, n_neighbors=10, n_pcs=20)
        toc = time.perf_counter()
        print(f"NBRS in {toc - tic:0.4f} seconds")
        #
        tic = time.perf_counter()
        sc.tl.leiden(axcat)
        toc = time.perf_counter()
        print(f"LEIDEN in {toc - tic:0.4f} seconds")
        #
        tic = time.perf_counter()
        sc.tl.umap(axcat)
        toc = time.perf_counter()
        print(f"UMAP in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    axcat.write_h5ad(axcat_file)
    toc = time.perf_counter()
    print(f"WRITE h5ad in {toc - tic:0.4f} seconds")
    #
    if plot_keys is not None:
        tic = time.perf_counter()
        sc.pl.umap(axcat, color=plot_keys,
                   save=axcat_plot_file)
        toc = time.perf_counter()
        print(f"PLOT in {toc - tic:0.4f} seconds")


def build_control_350k_ad(hd_dir, ngenes):
    adx = load_10X_180K_data()
    print(adx)
    _, cdx = load_covid_control_data(hd_dir, ngenes)
    print(cdx)
    lsx = [adx, cdx]
    axcat = prep_combat_intersect(lsx)
    return lsx, axcat


def run_10X_control_50K(gxtype, ngenes, pp_only, output_dir=".",
                        scement=True, post=False, workflow=False, filter_pct=0.0):
    axcat = load_10X_50K_data()
    print(axcat)
    infix_str = "10X_ctrl_scement_50k_"
    if not scement:
        infix_str = "sc_" + infix_str
    axcat_file = output_dir + "/" + infix_str + gxtype + ".h5ad"
    axcat_plot_file = "_" + infix_str + gxtype + ".png"
    pp_infix_str = "10X_ctrl_pp_50k_"
    pp_file = output_dir + "/" + pp_infix_str + gxtype + ".h5ad"
    if pp_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ntopg=ngenes,
                        scement=scement, plot_keys=None,
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_10X_control_350K(gxtype, ngenes, pp_only, output_dir=".",
                         scement=True, post=False, workflow=False, filter_pct=0.0):
    axcat = load_10X_350K_data()
    print(axcat)
    infix_str = "10X_ctrl_scement_350k_"
    if not scement:
        infix_str = "sc_" + infix_str
    axcat_file = output_dir + "/" + infix_str + gxtype + ".h5ad"
    axcat_plot_file = "_" + infix_str + gxtype + ".png"
    pp_infix_str = "10X_ctrl_pp_350k_"
    pp_file = output_dir + "/" + pp_infix_str + gxtype + ".h5ad"
    if pp_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=None, plot_keys=None,
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_control_350K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                     scement=True, post=False, workflow=False, filter_pct=0.0):
    _, axcat = build_control_350k_ad(hd_dir, ngenes)
    print(axcat)
    infix_str = "ctrl_scement_350k_"
    if not scement:
        infix_str = "sc_" + infix_str
    axcat_file = output_dir + "/" + infix_str + gxtype + ".h5ad"
    axcat_plot_file = "_" + infix_str + gxtype + ".png"
    if scement:
        axcat_file = output_dir + "/" + "ctrl_scement_350k_" + gxtype + ".h5ad"
        axcat_plot_file = "_ctrl_scement_350k_" + gxtype + ".png"
    else:
        axcat_file = output_dir + "/" + "sc_ctrl_combat_350k_" + gxtype + ".h5ad"
        axcat_plot_file = "_sc_ctrl_combat_350k_" + gxtype + ".png"
    pp_infix_str = "ctrl_pp_350k_"
    pp_file = output_dir + "/" + pp_infix_str + gxtype + ".h5ad"
    if pp_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, plot_keys=None, post=post,
                        workflow=workflow, filter_pct=filter_pct)


def build_infected_10k_ad(hd_dir, ngenes, union=False):
    ds10 = ['S-S062', 'S-M010-3', 'S-S054']
    print("Loading data from : ", ds10)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds10, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds10, ngenes)
    return lsx, axcat


def build_infected_25k_ad(hd_dir, ngenes, union=False):
    ds25 = ['S-S035-1', 'S-S037', 'S-S013', 'S-M077', 'S-S036-3']
    print("Loading data from : ", ds25)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds25, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds25, ngenes)
    return lsx, axcat


def build_infected_50k_ad(hd_dir, ngenes, union=False):
    ds50 = ['S-M018', 'S-S050', 'S-M066', 'S-S013', 'S-S035-1',
            'S-S001-2', 'S-S081', 'S-M037']
    print("Loading data from : ", ds50)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds50, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds50, ngenes)
    return lsx, axcat


def build_infected_150k_ad(hd_dir, ngenes, union=False):
    ds150 = ['S-M061-2', 'S-M041-1', 'S-S021-4', 'S-M044-1', 'S-M043-2',
             'S-M053', 'S-S018', 'S-S036-3', 'S-S052', 'S-S035-3',
             'S-M056', 'S-S087-2', 'S-M049', 'S-M020', 'S-M001',
             'S-S082', 'S-M035-1', 'S-M012', 'S-S083', 'S-S050',
             'S-S027', 'S-M018', 'S-S086-2', 'S-S061']
    print("Loading data from : ", ds150)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds150, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds150, ngenes)
    return lsx, axcat


def build_infected_350k_ad(hd_dir, ngenes, union=False):
    ds350 = [
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
    ]
    print("Loading data from : ", ds350)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds350, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds350, ngenes)
    return lsx, axcat


def build_infected_500k_ad(hd_dir, ngenes, union=False):
    ds500 = [
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
    ]
    print("Loading data from : ", ds500)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds500, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds500, ngenes)
    return lsx, axcat


def build_infected_700k_ad(hd_dir, ngenes, union=False):
    ds700 = [
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
    ]
    print("Loading data from : ", ds700)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds700, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds700, ngenes)
    return lsx, axcat


def build_infected_full_ad(hd_dir, ngenes, union=False):
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, None, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, None, ngenes)
    return lsx, axcat


def build_covid_full_ad(hd_dir, ngenes, union=False):
    if union:
        lsx, axcat = load_covid_full_data_union(hd_dir, None, ngenes)
    else:
        lsx, axcat = load_covid_full_data(hd_dir, None, ngenes)
    return lsx, axcat


def build_infected_50k_ad_union(hd_dir, ngenes):
    ds50 = ['S-M018', 'S-S050', 'S-M066', 'S-S013', 'S-S035-1',
            'S-S001-2', 'S-S081', 'S-M037']
    print("Loading data from : ", ds50)
    lsx, axcat = load_covid_infected_data_union(hd_dir, ds50, ngenes)
    return lsx, axcat


def get_sub_file_name(output_dir, cxtype, gxtype, scement,
                      union, dis_flag, filter_pct):
    if dis_flag:
        if scement:
            file_name = "dis_scement"
        else:
            file_name = "sc_dis_combat"
        pp_file_name = "dis_pp"
    else:
        if scement:
            file_name = "cvd_scement"
        else:
            file_name = "sc_cvd_combat"
        pp_file_name = "cvd_pp"
    file_name = file_name + "_" + cxtype + "_" + gxtype
    if dis_flag:
        pp_file_name = pp_file_name + "_" + cxtype + "_" + gxtype
    else:
        pp_file_name = pp_file_name + "_" + cxtype + "_" + gxtype
    if union:
        file_name = file_name + "_union"
        pp_file_name = pp_file_name + "_union"
    else:
        file_name = file_name + "_intersect"
        pp_file_name = pp_file_name + "_intersect"
    if filter_pct > 0:
        file_name = file_name + "_filter" + str(filter_pct * 100)
    axcat_file = output_dir + "/" + file_name + ".h5ad"
    pp_axcat_file = output_dir + "/" + pp_file_name + ".h5ad"
    axcat_plot_file = "_" + file_name + ".png"
    return axcat_file, axcat_plot_file, pp_axcat_file


def get_grow_file_name(output_dir, cxtype, gxtype, subxtype,
                       scement, union, dis_flag, filter_pct):
    if dis_flag:
        if scement:
            file_name = "dis_scement"
        else:
            file_name = "sc_dis_combat"
        pp_file_name = "dis_pp"
    else:
        if scement:
            file_name = "cvd_scement"
        else:
            file_name = "sc_cvd_combat"
        pp_file_name = "cvd_pp"
    file_name = file_name + "_" + cxtype + "_" + gxtype
    if dis_flag:
        pp_file_name = pp_file_name + "_" + cxtype + "_" + gxtype
    else:
        pp_file_name = pp_file_name + "_" + cxtype + "_" + gxtype
    file_name = file_name + "_g" + subxtype
    pp_file_name = pp_file_name + "_g" + subxtype
    if union:
        file_name = file_name + "_union"
        pp_file_name = pp_file_name + "_union"
    else:
        file_name = file_name + "_intersect"
        pp_file_name = pp_file_name + "_intersect"
    if filter_pct > 0:
        file_name = file_name + "_filter" + str(filter_pct * 100)
    axcat_file = output_dir + "/" + file_name + ".h5ad"
    pp_axcat_file = output_dir + "/" + pp_file_name + ".h5ad"
    axcat_plot_file = "_" + file_name + ".png"
    return axcat_file, axcat_plot_file, pp_axcat_file


def run_infected_10K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                     scement=True, post=False, union=False, workflow=False,
                     filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "10k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_10k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_25K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                     scement=True, post=False, union=False, workflow=False,
                     filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "25k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_25k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_50K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                     scement=True, post=False, union=False, workflow=False,
                     filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "50k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_50k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_150K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                      scement=True, post=False, union=False, workflow=False,
                      filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "150k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_150k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_350K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                      scement=True, post=False, union=False, workflow=False,
                      filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "350k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_350k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_500K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                      scement=True, post=False, union=False, workflow=False,
                      filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "500k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_500k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_700K(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                      scement=True, post=False, union=False, workflow=False,
                      filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "700k", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_700k_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_infected_full(hd_dir, gxtype, ngenes, pp_only, output_dir=".",
                      scement=True, post=False, union=False, workflow=False,
                      filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "full", gxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_infected_full_ad(hd_dir, ngenes, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:   # TODO:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def run_covid_full(hd_dir, gxtype, pp_only, output_dir=".",
                   scement=True, post=False, union=False, workflow=False,
                   filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_sub_file_name(
        output_dir, "cvd_full", gxtype, scement, union, False, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_covid_full_ad(hd_dir, 1000, union)
    print("AXCAT", axcat)
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:   # TODO:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        1000, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def main_covid_infected(hd_dir, cxtype, gxtype, ngenes, pp_only, output_dir=".",
                        scement=False, post=False,
                        union=False, workflow=False, filter_pct=0.0):
    if cxtype == "10K":
        run_infected_10K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                         scement, post, union, workflow, filter_pct)
    elif cxtype == "25K":
        run_infected_25K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                         scement, post, union, workflow, filter_pct)
    elif cxtype == "50K":
        run_infected_50K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                         scement, post, union, workflow, filter_pct)
    elif cxtype == "150K":
        run_infected_150K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                          scement, post, union, workflow, filter_pct)
    elif cxtype == "350K":
        run_infected_350K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                          scement, post, union, workflow, filter_pct)
    elif cxtype == "500K":
        run_infected_500K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                          scement, post, union, workflow, filter_pct)
    elif cxtype == "700K":
        run_infected_700K(hd_dir, gxtype, ngenes, pp_only, output_dir,
                          scement, post, union, workflow, filter_pct)
    elif cxtype == "full":
        run_infected_full(hd_dir, gxtype, ngenes, pp_only, output_dir,
                          scement, post, union, workflow, filter_pct)


def main_control_other(hd_dir, cxtype, gxtype, ngenes, pp_only, output_dir=".",
                       scement=False, post=False,
                       union=False, workflow=False, filter_pct=0.0):
    if cxtype == "50K":
        run_10X_control_50K(gxtype, ngenes, pp_only, output_dir,
                            scement, post, workflow, filter_pct)
    elif cxtype == "350K":
        # run_control_350K(hd_dir, gxtype, ngenes, output_dir, scement)
        run_10X_control_350K(gxtype, ngenes, pp_only, output_dir,
                             scement, post, workflow, filter_pct)
    elif cxtype == "full":
        run_covid_full(hd_dir, gxtype, pp_only, output_dir,
                       scement, post, union, workflow, filter_pct)
    else:
        print("Not implemented")


def main_sub(hd_dir, cxtype, gxtype, ngenes, pp_only, output_dir=".",
             scement=False, infected=False, post=False,
             union=False, workflow=False, filter_pct=0.0):
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    if infected:
        main_covid_infected(hd_dir, cxtype, gxtype, ngenes, pp_only, output_dir,
                            scement, post, union, workflow, filter_pct)
    else:
        main_control_other(hd_dir, cxtype, gxtype, ngenes, pp_only, output_dir,
                           scement, post, union, workflow, filter_pct)


def build_grw_150k_ad(hd_dir, subxtype, ngenes, union=False):
    ds150 = []
    if subxtype == "25K":
        ds150 = sd.subset150_25k()
    elif subxtype == "50K":
        ds150 = sd.subset150_50k()
    elif subxtype == "60K":
        ds150 = sd.subset150_60k()
    elif subxtype == "75K":
        ds150 = sd.subset150_75k()
    elif subxtype == "80K":
        ds150 = sd.subset150_80k()
    elif subxtype == "100K":
        ds150 = sd.subset150_100k()
    elif subxtype == "125K":
        ds150 = sd.subset150_125k()
    elif subxtype == "150K":
        ds150 = ['S-M061-2', 'S-M041-1', 'S-S021-4', 'S-M044-1', 'S-M043-2',
                 'S-M053', 'S-S018', 'S-S036-3', 'S-S052', 'S-S035-3',
                 'S-M056', 'S-S087-2', 'S-M049', 'S-M020', 'S-M001',
                 'S-S082', 'S-M035-1', 'S-M012', 'S-S083', 'S-S050',
                 'S-S027', 'S-M018', 'S-S086-2', 'S-S061']
    else:
        return None, None
    print("Loading data from : ", ds150)
    if union:
        lsx, axcat = load_covid_infected_data_union(hd_dir, ds150, ngenes)
    else:
        lsx, axcat = load_covid_infected_data(hd_dir, ds150, ngenes)
    return lsx, axcat


def run_grow_150K(hd_dir, gxtype, subxtype, ngenes, pp_only, output_dir=".",
                  scement=True, post=False, union=False, workflow=False,
                  filter_pct=0.0):
    axcat_file, axcat_plot_file, pp_axcat_file = get_grow_file_name(
        output_dir, "150K", gxtype, subxtype, scement, union, True, filter_pct)
    print("Output : ", axcat_file, axcat_plot_file)
    _, axcat = build_grw_150k_ad(hd_dir, subxtype, ngenes, union)
    print("AXCAT", axcat)
    if axcat is None:
        print(subxtype, " not supported for 150K")
        return None
    if pp_only:
        axcat.write_h5ad(pp_axcat_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file,
                        ngenes, scement, cvats=['batch', 'type'],
                        plot_keys=['pid', 'covid', 'cell_type'],
                        post=post, workflow=workflow, filter_pct=filter_pct)


def main_grow(hd_dir, cxtype, gxtype, subxtype, ngenes, pp_only, output_dir=".",
              scement=False, post=False, union=False,
              workflow=False, filter_pct=0.0):
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    if cxtype == "150K":
        run_grow_150K(hd_dir, gxtype, subxtype, ngenes, pp_only, output_dir,
                      scement, post, union, workflow, filter_pct)
    else:
        print("Not Implemented")


def ipython_info():
    ip = False
    if 'ipykernel' in sys.modules:
        ip = 'notebook'
    elif 'IPython' in sys.modules:
        ip = 'terminal'
    return ip


if __name__ == "__main__" and ipython_info() is False:
    rctypes = ["10K", "25K", "50K", "150K", "350K", "500K", "700K", "full"]
    default_in_dir = "./h5ad_full/"
    default_out_dir = "./"
    ngenes = {"500": 500, "1K": 1000, "5K": 5000, "10K": 10000, "full": None}
    rgtypes = list(ngenes.keys())
    #
    parser = argparse.ArgumentParser(prog="CVD",
                                     description="CVD COMBAT",
                                     epilog="CVD HELP")
    subparsers = parser.add_subparsers(help="Either Subset mode", dest="which")
    #
    #
    parser_sub = subparsers.add_parser("sub", help="Random subsets from full data")
    parser_sub.add_argument("cxtype", choices=rctypes)
    parser_sub.add_argument("gxtype", choices=rgtypes)
    parser_sub.add_argument("-d", "--infected", action="store_true")
    parser_sub.add_argument("-i", "--in_dir", default=default_in_dir)
    parser_sub.add_argument("-o", "--out_dir", default=default_out_dir)
    parser_sub.add_argument("-p", "--post", action="store_true")
    parser_sub.add_argument("-s", "--scement", action="store_true")
    parser_sub.add_argument("-u", "--union", action="store_true")
    parser_sub.add_argument("-w", "--workflow", action="store_true")
    parser_sub.add_argument("-y", "--pp_only", action="store_true")
    parser_sub.add_argument("-f", "--filter_pct", type=float, default=0.0)
    #
    #
    parser_grw = subparsers.add_parser("grw", help="Randomly grown datasets from full data")
    rctypes = ["150K"]
    parser_grw.add_argument("cxtype", choices=rctypes)
    parser_grw.add_argument("gxtype", choices=rgtypes)
    parser_grw.add_argument("subxtype", type=str)
    parser_grw.add_argument("-i", "--in_dir", default=default_in_dir)
    parser_grw.add_argument("-o", "--out_dir", default=default_out_dir)
    parser_grw.add_argument("-p", "--post", action="store_true")
    parser_grw.add_argument("-s", "--scement", action="store_true")
    parser_grw.add_argument("-u", "--union", action="store_true")
    parser_grw.add_argument("-w", "--workflow", action="store_true")
    parser_grw.add_argument("-y", "--pp_only", action="store_true")
    parser_grw.add_argument("-f", "--filter_pct", type=float, default=0.0)
    #
    #
    in_args = parser.parse_args(sys.argv[1:])
    print(in_args)
    if in_args.which == "sub":
        main_sub(in_args.in_dir, in_args.cxtype, in_args.gxtype,
                 ngenes[in_args.gxtype], in_args.pp_only, in_args.out_dir,
                 in_args.scement, in_args.infected, in_args.post,
                 in_args.union, in_args.workflow, in_args.filter_pct)
    if in_args.which == "grw":
        print("Grow sub ")
        main_grow(in_args.in_dir, in_args.cxtype, in_args.gxtype,
                  in_args.subxtype, ngenes[in_args.gxtype], in_args.pp_only,
                  in_args.out_dir, in_args.scement, in_args.post,
                  in_args.union, in_args.workflow, in_args.filter_pct)
