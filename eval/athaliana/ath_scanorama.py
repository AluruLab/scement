import pathlib
import sys
import time
import numpy as np
import anndata as an
import scanpy as sc
import scanorama as scanr
from pathlib import PurePath

DATA_DIR = "../../data/athaliana/"
ANALYSIS_OUT_DIR = "./athaliana_out/"


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


def scanorama_integrate2(scn_adatas, out_dir, file_sfx, plots=["batch", "cell_type"]):
    tic = time.perf_counter()
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    axcat_plot_file = "scanorama_" + file_sfx + ".png"
    axcat_results_file = out_dir + "/scanorama_" + file_sfx + ".h5ad"
    scn_adatas_cor = scanr.correct_scanpy(scn_adatas, return_dimred=True)
    scanr.integrate_scanpy(scn_adatas_cor, dimred=50)
    axcat = sc.concat(scn_adatas_cor, uns_merge="unique")
    sc.pp.neighbors(axcat, use_rep="X_scanorama")
    sc.tl.umap(axcat)
    sc.tl.leiden(axcat, key_added="clusters")
    sc.pl.umap(axcat, color=plots, save=axcat_plot_file)
    axcat.write(PurePath(axcat_results_file))
    toc = time.perf_counter()
    print(f"SCANORMA-INTEGRATE V2 in {toc - tic:0.4f} seconds")
    return axcat


def scanorama_integrate(
    axcat, out_dir, file_sfx, plots=["batch", "stress", "cell_type"]
):
    tic = time.perf_counter()
    sc._settings.ScanpyConfig.figdir = pathlib.Path(out_dir)
    axcat_plot_file = "scanorama1_" + file_sfx + ".png"
    axcat_results_file = out_dir + "/scanorama1_" + file_sfx + ".h5ad"
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
    sc.pl.umap(axcat, color=plots, save=axcat_plot_file)
    axcat.write(PurePath(axcat_results_file))
    toc = time.perf_counter()
    print(f"SCANORMA-INTEGRATE V1 in {toc - tic:0.4f} seconds")
    return axcat


def run_wt_gt_2dsets(data_dir, output_dir):
    tic = time.perf_counter()
    t61lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-158761/E-GEOD-158761-LB.h5ad")
    t19lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-121619/E-GEOD-121619-LB.h5ad")
    # 1619LB
    print("121619 LB : ", t19lblad.shape, type(t19lblad.X), end=" ")
    filter_scale(t19lblad)
    print(t19lblad.shape, type(t19lblad.X))
    # 8761LB
    print("158761 LB : ", t61lblad.shape, type(t61lblad.X), end=" ")
    filter_scale(t61lblad)
    print(t61lblad.shape, type(t61lblad.X))

    full_ads = [t19lblad, t61lblad]
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    toc = time.perf_counter()
    print(f"PREP TIME in {toc - tic:0.4f} seconds")
    # wt_axcat_norm_file = output_dir + "/" + "scanorama_norm_ath_wt2dsets.h5ad"
    # wt_ad_objects = [t19lblad, t61lblad]
    # wt_axcat = prep_intersect(wt_ad_objects)
    # wt_axcat.write(PurePath(wt_axcat_norm_file))
    # scanorama_integrate2(wt_axcat, output_dir, "ath_wt2dsets")
    scanorama_integrate2(
        full_ads, output_dir, "ath_wt2dsets", plots=["batch", "cell_type"]
    )


def run_copilot_dsets(data_dir, output_dir):
    tic = time.perf_counter()
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

    full_ads = [an.read_h5ad(data_dir + "/GSE152766/" + bx + ".h5ad") for bx in batches]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    toc = time.perf_counter()
    print(f"PREP TIME in {toc - tic:0.4f} seconds")
    #    cp_axcat_norm_file = output_dir + "/" + "scanorama_ath_coplit_norm.h5ad"
    #    cp_axcat = prep_intersect(full_ads)
    #    cp_axcat.write(PurePath(cp_axcat_norm_file))
    #    scanorama_integrate(cp_axcat, output_dir, "ath_copilot_",
    #                        plots=['batch', 'cell_type'])
    scanorama_integrate2(
        full_ads, output_dir, "ath_copilot", plots=["batch", "cell_type"]
    )


def run_filter1_jb_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761-FLT1.h5ad")
    jb_ds = an.read_h5ad(data_dir + "./E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT1.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    toc = time.perf_counter()
    print(f"PREP TIME in {toc - tic:0.4f} seconds")
    #    cp_axcat = prep_intersect(full_ads)
    #    cp_axcat_norm_file = output_dir + "/" + "scanorama_ath_jb_gala_flt1_norm.h5ad"
    #    cp_axcat.write(cp_axcat_norm_file)
    #    scanorama_integrate(cp_axcat, output_dir, "ath_jb_gala_flt_",
    #                        plots=['batch', 'stress', 'cell_type'])
    scanorama_integrate2(
        full_ads, output_dir, "ath_jb_gala_flt1", plots=["batch", "stress", "cell_type"]
    )


def run_filter2_jb_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761-FLT1.h5ad")
    jb_ds = an.read_h5ad(data_dir + "./E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT2.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    toc = time.perf_counter()
    print(f"PREP TIME in {toc - tic:0.4f} seconds")
    #    cp_axcat = prep_intersect(full_ads)
    #    cp_axcat_norm_file = output_dir + "/" + "scanorama_ath_jb_gala_flt1_norm.h5ad"
    #    cp_axcat.write(cp_axcat_norm_file)
    #    scanorama_integrate(cp_axcat, output_dir, "ath_jb_gala_flt_",
    #                        plots=['batch', 'stress', 'cell_type'])
    scanorama_integrate2(
        full_ads, output_dir, "ath_jb_gala_flt2", plots=["batch", "stress", "cell_type"]
    )


def run_jb_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761.h5ad")
    jb_ds = an.read_h5ad(data_dir + "./E-GEOD-121619/E-GEOD-121619-FULL-LB.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    toc = time.perf_counter()
    print(f"PREP TIME in {toc - tic:0.4f} seconds")
    #    cp_axcat = prep_intersect(full_ads)
    #    cp_axcat_norm_file = output_dir + "/" + "scanorama_ath_jb_gala_norm.h5ad"
    #    cp_axcat.write(cp_axcat_norm_file)
    #    scanorama_integrate(cp_axcat, output_dir, "ath_jb_gala_",
    #                        plots=['batch', 'cell_type'])
    scanorama_integrate2(
        full_ads, output_dir, "ath_jb_gala", plots=["batch", "cell_type"]
    )


def run_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761.h5ad")
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    #    filter_scale(gala_ds)
    #    cp_axcat = prep_intersect(gala_ds)
    #    cp_axcat_norm_file = output_dir + "/" + "scanorama_ath_jb_gala_norm.h5ad"
    #    cp_axcat.write(cp_axcat_norm_file)
    #    scanorama_integrate(cp_axcat, output_dir, "ath_gala_",
    #                        plots=['batch', 'cell_type'])
    batches = list(set(gala_ds.obs["batch"].tolist()))
    print(batches)
    adatas = []
    for batch in batches:
        adatas.append(gala_ds[gala_ds.obs["batch"] == batch, ])  # type:ignore
    #
    toc = time.perf_counter()
    print(f"PREP TIME in {toc - tic:0.4f} seconds")
    scanorama_integrate2(adatas, output_dir, "ath_gala", plots=["batch", "cell_type"])


def main(rtype):
    if rtype == "2ds":
        run_wt_gt_2dsets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "copilot":
        run_copilot_dsets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "jb_gala":
        run_jb_gala_datasets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "gala":
        run_gala_datasets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "jb_gala_flt1":
        run_filter1_jb_gala_datasets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "jb_gala_flt2":
        run_filter2_jb_gala_datasets(DATA_DIR, ANALYSIS_OUT_DIR)
    else:
        run_filter2_jb_gala_datasets(DATA_DIR, ANALYSIS_OUT_DIR)


if __name__ == "__main__":
    run_choices = ["2ds", "copilot", "jb_gala", "gala", "jb_gala_flt1",
                   "jb_gala_flt2", "default"]
    run_option = "/".join(run_choices)
    if len(sys.argv) <= 1:
        main("default")
    else:
        rtype = sys.argv[1]
        main(rtype)
