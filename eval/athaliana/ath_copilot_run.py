import pathlib
import time
import numpy as np

# import pandas as pd
import anndata as an
import scanpy as sc
import scanorama as scanr

DATA_DIR = "./ath_integ/"
OUTPUT_DIR = "./athaliana/"


def filter_scale(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)


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


def scanorama_integrate(
    axcat, out_dir, file_pfx, plots=["batch", "stress", "cell_type"]
):
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    axcat_plot_file = file_pfx + "_scanorama_umap.png"
    axcat_results_file = out_dir + "/" + file_pfx + "_scanorama.h5ad"
    # split per batch into new objects.
    print(axcat.obs["batch"])
    batches = list(set(axcat.obs["batch"].tolist()))
    print(batches)
    adatas_map = {}
    for batch in batches:
        adatas_map[batch] = axcat[axcat.obs["batch"] == batch,]
    scn_adatas = adatas_map.values()
    scanr.integrate_scanpy(scn_adatas, dimred=50)
    # Get all the integrated matrices.
    scanorama_int = [ad.obsm["X_scanorama"] for ad in scn_adatas]
    # make into one matrix.
    all_s = np.concatenate(scanorama_int)
    print(all_s.shape)
    # add to the AnnData object
    axcat.obsm["Scanorama"] = all_s
    # tsne and umap
    sc.pp.neighbors(axcat, n_pcs=50, use_rep="Scanorama")
    # sc.tl.umap(adata)
    # sc.tl.tsne(adata, n_pcs = 50, use_rep = "Scanorama")
    sc.tl.leiden(axcat)
    sc.tl.umap(axcat)
    sc.pl.umap(
        axcat,
        color=plots,  # ['batch', 'cell_type', 'cluster_ext_type'],
        save=axcat_plot_file,
    )
    axcat.write(axcat_results_file)
    return axcat


def main():
    gala_ds = an.read_h5ad(DATA_DIR + "/GSE158761/GSE158761.h5ad")
    jb_ds = an.read_h5ad(DATA_DIR + "./E-GEOD-121619/E-GEOD-121619-FULL-LB.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(OUTPUT_DIR)
    cp_axcat = prep_intersect(full_ads)
    cp_axcat_norm_file = OUTPUT_DIR + "/scanorama_ath_jb_gala_norm.h5ad"
    cp_axcat.write(pathlib.PurePath(cp_axcat_norm_file))
    scanorama_integrate(cp_axcat, OUTPUT_DIR, "ath_copilot_", plots=["batch", "cell_type"])
    #
    batches = [
        "sc_1",
        "sc_31",
        "tnw2",
        "sc_11",
        "sc_51",
        "sc_37",
        "sc_9_at",
        "tnw1",
        "sc_40",
        "sc_12",
        "col0",
        "sc_30",
        "sc_10_at",
    ]
    #
    full_ads = [an.read_h5ad(DATA_DIR + "/GSE152766/" + bx + ".h5ad") for bx in batches]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(OUTPUT_DIR)
    cp_axcat_norm_file = OUTPUT_DIR + "/scanorama_ath_coplit_norm.h5ad"
    cp_axcat = prep_intersect(full_ads)
    cp_axcat.write(pathlib.PurePath(cp_axcat_norm_file))
    scanorama_integrate(cp_axcat, OUTPUT_DIR, "ath_copilot_", plots=["batch", "cell_type"])


if __name__ == "__main__":
    main()
