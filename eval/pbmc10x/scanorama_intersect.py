import anndata as ad
import scanpy as sc
import scanorama as scanr
import argparse
import pathlib
import json
import time


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
    axcat.write(pathlib.PurePath(axcat_results_file))
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


def main_json(json_data, out_dir, out_file_pfx):
    data_dir = json_data['DATADIR']
    in_files = json_data['ADFILE']
    # batch_names = json_data['BATCH']
    # out_dir = json_data['OUTDIR']
    # out_file = json_data['OUTFILE']
    # out_file_pfx = pathlib.Path(out_file).stem
    # results_file = out_dir + "/" + out_file_pfx + ".h5ad"
    #
    tic = time.perf_counter()
    ad_objects = [ad.read_h5ad(data_dir + "/" + fx) for fx in in_files]
    for x in ad_objects:
        x.var.index = [y for y in x.var.gene_ids]
    au_objects = [ax[:, ~ax.var.gene_ids.duplicated()] for ax in ad_objects]
    toc = time.perf_counter()
    print(f"LOAD in {toc - tic:0.4f} seconds")
    print(ad_objects)
    print(au_objects)
    #
    #
    tic = time.perf_counter()
    adatas = au_objects
    tx = scanorama_integrate2(adatas, out_dir, out_file_pfx,
                              None, None, False)
    # tx = scanorama.correct_scanpy(adatas, return_dimred=True)
    toc = time.perf_counter()
    print(f"CORRECT in {toc - tic:0.4f} seconds")
    print(tx)


if __name__ == "__main__":
    help_str = """Path to input json file, File expected to have following format: """
    parser = argparse.ArgumentParser(description='Integrate w. scanorama.')
    parser.add_argument('json_file', type=str, help=help_str)
    parser.add_argument("out_dir")
    parser.add_argument("out_file_prefix")
    args = parser.parse_args()
    with open(args.json_file) as f:
        json_data = json.load(f)
    #
    main_json(json_data, args.out_dir, args.out_file_prefix)
