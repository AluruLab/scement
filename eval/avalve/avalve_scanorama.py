import pathlib
import time
import numpy as np
import anndata as an
import scanpy as sc
import scanorama as scanr
from pathlib import PurePath

# Generate analyses with Scanorama
# source data is from DATA_DIR, and output is ANALYSIS_RESULTS_DIR
DATA_DIR = "../../data/aortic_valve/"
ANALYSIS_RESULTS_DIR = "./avalve_out/"


def filter_scale(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)
    # sc.pp.scale(adata, max_value=10)


def prep_intersect(ad_objects):
    tic = time.perf_counter()
    au_objects = [ax[:, ~ax.var.gene_ids.duplicated()] for ax in ad_objects]
    print(au_objects)
    #
    # axcat = ad.concat(au_objects)
    axcat = an.concat(au_objects, merge="same")
    axcat.obs_names_make_unique()
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    return axcat


def save_axcat(axcat, out_dir, file_sfx, plots):
    tic = time.perf_counter()
    axcat_plot_file = "scanorama_" + file_sfx + ".png"
    axcat_results_file = out_dir + "/scanorama_" + file_sfx + ".h5ad"
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    sc.pl.umap(axcat, color=plots, save=axcat_plot_file)
    axcat.write(PurePath(axcat_results_file))
    axcat = axcat[axcat.obs["cell_type"] != "Unknown"]
    filter_axcat_plot_file = "scanorama_" + file_sfx + "_filter.png"
    sc.pl.umap(
        axcat,
        color=plots,  # ['batch', 'cell_type', 'cluster_ext_type'],
        save=filter_axcat_plot_file,
    )
    sfilter_axcat_plot_file = "scanorama_" + file_sfx + "_sfilter.png"
    sc.pl.umap(axcat, color=["batch", "cell_stype"], save=sfilter_axcat_plot_file)
    axcat2 = axcat[axcat.obs["cell_type"] != "Monocytes"]
    sfilter2_axcat_plot_file = "scanorama_" + file_sfx + "_final.png"
    sc.pl.umap(axcat2, color=["batch", "cell_stype"], save=sfilter2_axcat_plot_file)
    toc = time.perf_counter()
    print(f"SAVE-AXCAT in {toc - tic:0.4f} seconds")
    return axcat2


def scanorama_integrate2(scn_adatas, out_dir, file_pfx, plots=["batch", "cell_type"]):
    tic = time.perf_counter()
    print(scn_adatas)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    scn_adatas_cor = scanr.correct_scanpy(scn_adatas, return_dimred=True)
    scanr.integrate_scanpy(scn_adatas_cor, dimred=50)
    axcat = sc.concat(scn_adatas_cor, uns_merge="unique")
    sc.pp.neighbors(axcat, use_rep="X_scanorama")
    sc.tl.umap(axcat)
    sc.tl.leiden(axcat, key_added="clusters")
    save_axcat(axcat, out_dir, file_pfx, plots)
    toc = time.perf_counter()
    print(f"SCANORMA-INTEGRATE V2 (w. save) in {toc - tic:0.4f} seconds")
    return axcat


def scanorama_integrate(axcat, out_dir, file_pfx, plots=["batch", "cell_type"]):
    tic = time.perf_counter()
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    # split per batch into new objects.
    print(axcat.obs["batch"])
    batches = list(set(axcat.obs["batch"].tolist()))
    print(batches)
    adatas_map = {}
    for batch in batches:
        adatas_map[batch] = axcat[axcat.obs["batch"] == batch, ]
    scn_adatas = adatas_map.values()
    scanr.integrate_scanpy(scn_adatas, dimred=50)
    # Get all the integrated matrices.
    scanorama_int = [ad.obsm["X_scanorama"] for ad in scn_adatas]
    # make into one matrix.
    all_s = np.concatenate(scanorama_int)
    print(all_s.shape)
    # add to the AnnData object
    axcat.obsm["X_scanorama"] = all_s
    # tsne and umap
    sc.pp.neighbors(axcat, n_pcs=50, use_rep="X_scanorama")
    # sc.tl.umap(adata)
    # sc.tl.tsne(adata, n_pcs = 50, use_rep = "Scanorama")
    sc.tl.leiden(axcat)
    sc.tl.umap(axcat)
    save_axcat(axcat, out_dir, file_pfx, plots)
    toc = time.perf_counter()
    print(f"SCANORMA-INTEGRATE V1 (w. save) in {toc - tic:0.4f} seconds")
    return axcat


def run_avalve_dsets1(data_dir, output_dir):
    tic = time.perf_counter()
    batches = ["H1", "H2", "C3", "C4"]
    condition = ["Healthy", "Healthy", "Calcified", "Calcified"]
    #
    full_ads = [an.read_h5ad(data_dir + "/" + bx + ".h5ad") for bx in batches]
    for ix in range(len(full_ads)):
        bx = batches[ix]
        full_ads[ix].obs["health"] = bx[0]
        full_ads[ix].obs["condition"] = condition[ix]
        full_ads[ix].obs["cell_stype"] = full_ads[ix].obs[
            "health"
        ].astype(  # type:ignore
            "str"
        ) + full_ads[
            ix
        ].obs[
            "cell_type"
        ].astype(
            "str"
        )  # type:ignore
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_norm_file = output_dir + "/" + "scanorama_normalized_avalve.h5ad"
    cp_axcat = prep_intersect(full_ads)
    cp_axcat.write(PurePath(cp_axcat_norm_file))
    toc = time.perf_counter()
    print(f"PREP DATA in {toc - tic:0.4f} seconds")
    scanorama_integrate(
        cp_axcat, output_dir, "v1_avalve", plots=["batch", "condition", "cell_stype"]
    )


def run_avalve_dsets2(data_dir, output_dir):
    tic = time.perf_counter()
    batches = ["H1", "H2", "C3", "C4"]
    condition = ["Healthy", "Healthy", "Calcified", "Calcified"]

    full_ads = [an.read_h5ad(data_dir + "/" + bx + ".h5ad") for bx in batches]
    for ix in range(len(full_ads)):
        bx = batches[ix]
        full_ads[ix].obs["health"] = bx[0]
        full_ads[ix].obs["condition"] = condition[ix]
        full_ads[ix].obs["cell_stype"] = full_ads[ix].obs[
            "health"
        ].astype(  # type:ignore
            "str"
        ) + full_ads[
            ix
        ].obs[
            "cell_type"
        ].astype(
            "str"
        )  # type:ignore
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_norm_file = output_dir + "/" + "scanorama_normalized_avalve.h5ad"
    cp_axcat = prep_intersect(full_ads)
    cp_axcat.write(PurePath(cp_axcat_norm_file))
    toc = time.perf_counter()
    print(f"PREP DATA in {toc - tic:0.4f} seconds")
    scanorama_integrate2(full_ads, output_dir, "avalve", plots=["batch", "cell_type"])


def main():
    run_avalve_dsets2(DATA_DIR, ANALYSIS_RESULTS_DIR)


if __name__ == "__main__":
    main()
