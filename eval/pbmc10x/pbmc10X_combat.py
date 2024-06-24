import pandas as pd
import scanpy as sc
import anndata as an
import numpy as np
import sys
import time
import pathlib
import argparse
import json
import scement as sct


def combat_and_plot(axcat, axcat_plot_file, axcat_file,
                    ntopg=5000, scement=True,
                    batch_key="batch",
                    # cvats=['batch', 'source'],
                    cvats=[],
                    plot_keys=['pid', 'covid', 'cell_type'],
                    run_dtype=np.float32,
                    filter_pct=0.0):
    if (axcat is not None) and (filter_pct > 0.0):
        tic = time.perf_counter()
        pre_batches = set(axcat.obs[batch_key])
        celldx, genedx = sc.pp.calculate_qc_metrics(axcat)
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
    if ntopg is not None:
        tic = time.perf_counter()
        sc.pp.highly_variable_genes(axcat, flavor="seurat", n_top_genes=ntopg)
        axcat = axcat[:, axcat.var.highly_variable]
        toc = time.perf_counter()
        print(f"HVG in {toc - tic:0.4f} seconds")
    print(axcat)
    print("Run Args: ", batch_key, cvats, run_dtype)
    tic = time.perf_counter()
    if scement:
        sct.sct_sparse(axcat, key=batch_key, covariates=cvats, inplace=True,
                       run_dtype=run_dtype)
    else:
        sct.combat(axcat, key=batch_key, covariates=cvats, inplace=True,
                   run_dtype=run_dtype)
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")
    #
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
    axcat.write_h5ad(axcat_file)
    #
    if plot_keys is not None:
        tic = time.perf_counter()
        sc.pl.umap(axcat, color=plot_keys,
                   save=axcat_plot_file)
        toc = time.perf_counter()
        print(f"PLOT in {toc - tic:0.4f} seconds")


def update_obs_var(adx,
                   gm_file="/project/spc/i3/eval/covid_atlas/meta/gene_mapping.csv"):
    new_obs = adx.obs
    new_obs['pid'] = new_obs['batch']
    adx.obs = new_obs
    lkdf = pd.read_csv(gm_file, index_col=0)
    lkdf2 = lkdf.copy()  # type: ignore
    lkdf2.index = lkdf2.gene_id
    adx.var.index = [y for y in adx.var.gene_ids]
    new_var = lkdf2.loc[adx.var.index,]
    new_var.index = new_var.gene_name
    new_var['gene_id2'] = new_var.gene_id
    new_var['gene_id'] = new_var.gene_name
    adx.var = new_var
    return adx


def load_10X_ad(in_ad_file,
                gm_file="/project/spc/i3/eval/covid_atlas/meta/gene_mapping.csv",
                no_update_var=False):
    tic = time.perf_counter()
    adx = an.read_h5ad(in_ad_file)
    if not no_update_var:
        adx = update_obs_var(adx, gm_file)
        adx.var_names_make_unique()
    toc = time.perf_counter()
    print(f"COMBO LOAD in {toc - tic:0.4f} seconds")
    return adx


def main_combo(in_ad_file, gxtype, ngenes, output_dir=".",
               scement=False, dtype=np.float32, no_update_var=False,
               filter_pct=0.0):
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    axcat = load_10X_ad(in_ad_file, no_update_var=no_update_var)
    print(axcat)
    filter_pfx = ""
    if filter_pct > 0:
        filter_pfx = "_filter" + str(int(filter_pct * 100))
    if scement:
        file_bname = "scement_" + in_ad_file.replace(".h5ad", "_") + gxtype + filter_pfx
    else:
        file_bname = "combat_" + in_ad_file.replace(".h5ad", "_") + gxtype + filter_pfx
    axcat_file = output_dir + "/" + file_bname + ".h5ad"
    axcat_plot_file = file_bname + ".png"
    combat_and_plot(axcat, axcat_plot_file, axcat_file,
                    ngenes, scement, cvats=None, plot_keys=None,
                    run_dtype=dtype, filter_pct=filter_pct)


def prep_intersect(json_data):
    data_dir = json_data['DATADIR']
    in_files = json_data['ADFILE']
    tic = time.perf_counter()
    ad_objects = [an.read_h5ad(data_dir + "/" + fx) for fx in in_files]
    for x in ad_objects:
        x.var.index = [y for y in x.var.gene_ids]
    ad_objects = [ax[:, ~ax.var.gene_ids.duplicated()] for ax in ad_objects]
    toc = time.perf_counter()
    print(f"LOAD in {toc - tic:0.4f} seconds")
    #
    tic = time.perf_counter()
    # axcat = ad.concat(au_objects)
    axcat = an.concat(ad_objects, merge='same')
    axcat.obs_names_make_unique()
    axcat = update_obs_var(axcat)
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    return axcat


def main_json(in_json, gxtype, ngenes, output_dir=".",
              scement=False, dtype=np.float32, filter_pct=0.0):
    with open(in_json) as f:
        json_data = json.load(f)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    axcat = prep_intersect(json_data)
    out_dir = json_data['OUTDIR']
    out_file = json_data['OUTFILE']
    axcat_file = out_dir + "/" + out_file.replace(".h5ad", "_" + gxtype + ".h5ad")
    axcat_plot_file = out_file.replace("h5ad", "png")
    combat_and_plot(axcat, axcat_plot_file, axcat_file,
                    ngenes, scement, cvats=None, plot_keys=None,
                    run_dtype=dtype, filter_pct=filter_pct)


if __name__ == "__main__":
    ngenes = {"500": 500, "1K": 1000, "5K": 5000, "10K": 10000, "full": None}
    dtypes = {"32": np.float32, "64": np.float64}
    rgtypes = list(ngenes.keys())
    parser = argparse.ArgumentParser(prog="CVD",
                                     description="CVD COMBAT",
                                     epilog="CVD HELP")
    parser.add_argument("in_combo")
    parser.add_argument("out_dir")
    parser.add_argument("gxtype", choices=rgtypes)
    parser.add_argument("-d", "--dtype", choices=list(dtypes.keys()), default="32")
    parser.add_argument("-s", "--scement", action="store_true")
    parser.add_argument("-n", "--no_update_var", action="store_true")
    parser.add_argument("-f", "--filter_pct", type=float, default=0.0)
    in_args = parser.parse_args(sys.argv[1:])
    print(in_args)
    if (in_args.in_combo.endswith("h5ad")):
        main_combo(in_args.in_combo, in_args.gxtype, ngenes[in_args.gxtype],
                   in_args.out_dir, in_args.scement, dtypes[in_args.dtype],
                   in_args.no_update_var, in_args.filter_pct)
    elif (in_args.in_combo.endswith("json")):
        main_json(in_args.in_combo, in_args.gxtype, ngenes[in_args.gxtype],
                  in_args.out_dir, in_args.scement, dtypes[in_args.dtype],
                  in_args.filter_pct)
