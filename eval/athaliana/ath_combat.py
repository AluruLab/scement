import scanpy as sc
import anndata as an
import sys
import time
import scement
import pathlib

DATA_DIR = "./ath_integ/"
ANALYSIS_OUT_DIR = "./athaliana/"


def filter_scale(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)


def prep_combat_intersect(ad_objects):
    tic = time.perf_counter()
    au_objects = [ax[:, ~ax.var.gene_ids.duplicated()] for ax in ad_objects]
    print(au_objects)
    axcat = an.concat(au_objects, merge='same')
    axcat.obs_names_make_unique()
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    return axcat


def run_wt_gt_2dsets(data_dir, output_dir):
    tic = time.perf_counter()
    t61lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-158761/E-GEOD-158761-LB.h5ad")
    t19lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-121619/E-GEOD-121619-LB.h5ad")
    wt_raw_axcat = prep_combat_intersect([t19lblad, t61lblad])
    wt_raw_axcat.write(output_dir + "scanpy_raw_ath_wt2dsets.h5ad")
    # 1619LB
    print("121619 LB : ", t19lblad.shape, type(t19lblad.X), end=" ")
    filter_scale(t19lblad)
    print(t19lblad.shape, type(t19lblad.X))
    # 8761LB
    print("158761 LB : ", t61lblad.shape, type(t61lblad.X), end=" ")
    filter_scale(t61lblad)
    print(t61lblad.shape, type(t61lblad.X))

    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    wt_axcat_plot_file = "combat_wt2dsets.png"
    wt_axcat_results_file = output_dir + "/" + "combat_wt2dsets.h5ad"
    wt_axcat_norm_file = output_dir + "/" + "scanpy_normalized_wt2dsets.h5ad"
    wt_ad_objects = [t19lblad, t61lblad]
    wt_axcat = prep_combat_intersect(wt_ad_objects)
    wt_axcat.write(wt_axcat_norm_file)
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    scement.combat(wt_axcat, key='batch', covariates=['age', 'stress'], inplace=True)
    sc.tl.pca(wt_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(wt_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(wt_axcat)
    sc.tl.umap(wt_axcat)
    sc.pl.umap(wt_axcat, color=['batch', 'stress', 'cell_type'],
               save=wt_axcat_plot_file)
    wt_axcat.write(wt_axcat_results_file)
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def run_wt_gt_4dsets(data_dir, output_dir):
    # t66wtad = an.read_h5ad("E-GEOD-152766/E-GEOD-152766-WT.h5ad")
    # t66btad = an.read_h5ad("E-GEOD-152766/E-GEOD-152766-BT.h5ad")
    # Raw combo files
    t19lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-121619/E-GEOD-121619-LB.h5ad")
    t13lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-123013/E-GEOD-123013-LB.h5ad")
    t13rhdad = an.read_h5ad(data_dir + "/" + "E-GEOD-123013/E-GEOD-123013-RHD.h5ad")
    t61lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-158761/E-GEOD-158761-LB.h5ad")
    wt_raw_axcat = prep_combat_intersect([t19lblad, t13lblad, t61lblad])
    wt_raw_axcat.write(output_dir + "scanpy_ath_wt_raw.h5ad")
    gt_raw_axcat = prep_combat_intersect([t19lblad, t13lblad, t13rhdad])
    gt_raw_axcat.write(output_dir + "scanpy_ath_gt_raw.h5ad")
    #
    t19lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-121619/E-GEOD-121619-LB.h5ad")
    print("121619 LB : ", t19lblad.shape, type(t19lblad.X), t19lblad.X.nnz, end=" ")  # type:ignore
    filter_scale(t19lblad)
    print(t19lblad.shape, type(t19lblad.X), t19lblad.X.nnz)  # type:ignore
    t19lblad.write(data_dir + "/" + "E-GEOD-121619/E-GEOD-121619-LB-NORM.h5ad")

    #
    t13lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-123013/E-GEOD-123013-LB.h5ad")
    print("123013 LB : ", t13lblad.shape, type(t13lblad), t13lblad.X.nnz, end=" ")  # type:ignore
    filter_scale(t13lblad)
    print(t13lblad.shape, type(t13lblad.X), t13lblad.X.nnz)  # type:ignore
    t13lblad.write(data_dir + "/" + "E-GEOD-123013/E-GEOD-123013-LB-NORM.h5ad")

    #
    t13rhdad = an.read_h5ad(data_dir + "/" + "E-GEOD-123013/E-GEOD-123013-RHD.h5ad")
    print("123013 RHD : ", t13rhdad.shape, type(t13rhdad), t13rhdad.X.nnz, end=" ")  # type:ignore
    filter_scale(t13rhdad)
    print(t13rhdad.shape, type(t13rhdad.X), t13rhdad.X.nnz)  # type:ignore
    t13rhdad.write(data_dir + "/" + "E-GEOD-123013/E-GEOD-123013-RHD-NORM.h5ad")
    #
    t61lblad = an.read_h5ad(data_dir + "/" + "E-GEOD-158761/E-GEOD-158761-LB.h5ad")
    print("158761 LB : ", t61lblad.shape, type(t61lblad), t61lblad.X.nnz, end=" ")  # type:ignore
    filter_scale(t61lblad)
    print(t61lblad.shape, type(t61lblad), t61lblad.X.nzz)  # type:ignore
    t61lblad.write(data_dir + "/" + "E-GEOD-158761/E-GEOD-158761-LB-NORM.h5ad")

    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    wt_axcat_plot_file = "combat_ath_wt.png"
    wt_axcat_results_file = output_dir + "/" + "combat_ath_wt.h5ad"
    wt_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_wt.h5ad"
    wt_ad_objects = [t19lblad, t13lblad, t61lblad]
    wt_axcat = prep_combat_intersect(wt_ad_objects)
    wt_axcat.write(wt_axcat_norm_file)
    scement.combat(wt_axcat, key='batch', covariates=['age', 'stress'], inplace=True)
    sc.tl.pca(wt_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(wt_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(wt_axcat)
    sc.tl.umap(wt_axcat)
    sc.pl.umap(wt_axcat, color=['batch', 'stress', 'cell_type'],
               save=wt_axcat_plot_file)
    wt_axcat.write(wt_axcat_results_file)

    gt_gxcat_plot_file = "combat_ath_wt.png"
    gt_gxcat_results_file = output_dir + "/" + "combat_ath_gt.h5ad"
    gt_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_gt.h5ad"
    gt_ad_objects = [t19lblad, t13lblad, t13rhdad]
    gt_axcat = prep_combat_intersect(gt_ad_objects)
    gt_axcat.write(gt_axcat_norm_file)
    scement.combat(gt_axcat, key='batch', covariates=['genotype', 'age', 'stress'], inplace=True)
    sc.tl.pca(gt_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(gt_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(gt_axcat)
    sc.tl.umap(gt_axcat)
    sc.pl.umap(gt_axcat, color=['batch', 'stress', 'cell_type'],
               save=gt_gxcat_plot_file)
    gt_axcat.write(gt_gxcat_results_file)


def run_col0_dsets(data_dir, output_dir):
    tic = time.perf_counter()
    t66col0ad = an.read_h5ad(data_dir + "/" + "E-GEOD-152766/E-GEOD-152766-COL0.h5ad")
    col0_axcat = t66col0ad
    print("152766 LB : ", col0_axcat.shape, type(col0_axcat), col0_axcat.X.nnz, end=" ")  # type:ignore
    filter_scale(col0_axcat)
    print(col0_axcat.shape, type(col0_axcat), col0_axcat.X.nnz)  # type:ignore
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    col0_axcat_plot_file = "combat_ath_66col0.png"
    # col0_axcat_results_file = output_dir + "/" + "combat_ath_66col0.h5ad"
    col0_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_66col0.h5ad"
    col0_axcat = col0_axcat[:, ~col0_axcat.var.gene_ids.duplicated()]
    col0_axcat.obs_names_make_unique()
    col0_axcat.write(col0_axcat_norm_file)
    print(col0_axcat)
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    #
    # axcat = ad.concat(au_objects)
    scement.combat(col0_axcat, key='batch', covariates=['genotype', 'age', 'stress'], inplace=True)
    sc.tl.pca(col0_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(col0_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(col0_axcat)
    sc.tl.umap(col0_axcat)
    sc.pl.umap(col0_axcat, color=['batch', 'stress', 'cell_type'], save=col0_axcat_plot_file)
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def run_copilot_dsets(data_dir, output_dir):
    batches = ['sc_1', 'sc_31', 'tnw2', 'sc_11', 'sc_51', 'sc_37',
               'sc_9_at', 'tnw1', 'sc_40', 'sc_12', 'col0', 'sc_30', 'sc_10_at']
    #
    tic = time.perf_counter()
    full_ads = [an.read_h5ad(data_dir + "/GSE152766/" + bx + ".h5ad") for bx in batches]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    cp_axcat = prep_combat_intersect(full_ads)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_plot_file = "combat_ath_copilot.png"
    cp_axcat_results_file = output_dir + "/" + "combat_ath_copilot.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_coplit.h5ad"
    cp_axcat = cp_axcat[:, ~cp_axcat.var.gene_ids.duplicated()]
    cp_axcat.obs_names_make_unique()
    cp_axcat.write(cp_axcat_norm_file)
    #
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    #
    scement.combat(cp_axcat, key='batch', inplace=True)
    sc.tl.pca(cp_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(cp_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(cp_axcat)
    sc.tl.umap(cp_axcat)
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_plot_file)
    cp_axcat.write(cp_axcat_results_file)
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def run_filter1_jb_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761-FLT1.h5ad")
    jb_ds = an.read_h5ad(data_dir + "./E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT1.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    cp_axcat = prep_combat_intersect(full_ads)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_plot_file = "combat_ath_jb_gala_flt1.png"
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala_flt1.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt1.h5ad"
    cp_axcat = cp_axcat[:, ~cp_axcat.var.gene_ids.duplicated()]
    cp_axcat.obs_names_make_unique()
    cp_axcat.write(pathlib.PurePath(cp_axcat_norm_file))
    #
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    #
    # axcat = ad.concat(au_objects)
    scement.combat(cp_axcat, key='batch', covariates=['age', 'stress'], inplace=True)
    sc.tl.pca(cp_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(cp_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(cp_axcat)
    sc.tl.umap(cp_axcat)
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_plot_file)
    cp_axcat.write(pathlib.PurePath(cp_axcat_results_file))
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def run_filter2_jb_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761-FLT1.h5ad")
    jb_ds = an.read_h5ad(data_dir + "./E-GEOD-121619/E-GEOD-121619-FULL-LB-FLT2.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    cp_axcat = prep_combat_intersect(full_ads)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_plot_file = "combat_ath_jb_gala_flt2.png"
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala_flt2.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    cp_axcat = cp_axcat[:, ~cp_axcat.var.gene_ids.duplicated()]
    cp_axcat.obs_names_make_unique()
    cp_axcat.write(pathlib.PurePath(cp_axcat_norm_file))
    #
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    #
    # axcat = ad.concat(au_objects)
    scement.combat(cp_axcat, key='batch', covariates=['age', 'stress'], inplace=True)
    sc.tl.pca(cp_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(cp_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(cp_axcat)
    sc.tl.umap(cp_axcat)
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_plot_file)
    cp_axcat.write(pathlib.PurePath(cp_axcat_results_file))
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def run_jb_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761.h5ad")
    jb_ds = an.read_h5ad(data_dir + "./E-GEOD-121619/E-GEOD-121619-FULL-LB.h5ad")
    full_ads = [gala_ds, jb_ds]
    for ix in range(len(full_ads)):
        filter_scale(full_ads[ix])
    cp_axcat = prep_combat_intersect(full_ads)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_plot_file = "combat_ath_jb_gala.png"
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala.h5ad"
    cp_axcat = cp_axcat[:, ~cp_axcat.var.gene_ids.duplicated()]
    cp_axcat.obs_names_make_unique()
    cp_axcat.write(pathlib.PurePath(cp_axcat_norm_file))
    #
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    #
    scement.combat(cp_axcat, key='batch', covariates=['age', 'stress'], inplace=True)
    sc.tl.pca(cp_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(cp_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(cp_axcat)
    sc.tl.umap(cp_axcat)
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_plot_file)
    cp_axcat.write(pathlib.PurePath(cp_axcat_results_file))
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def run_gala_datasets(data_dir, output_dir):
    tic = time.perf_counter()
    gala_ds = an.read_h5ad(data_dir + "/GSE158761/GSE158761.h5ad")
    filter_scale(gala_ds)
    cp_axcat = prep_combat_intersect(gala_ds)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_plot_file = "combat_ath_gala.png"
    cp_axcat_results_file = output_dir + "/" + "combat_ath_gala.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_gala.h5ad"
    cp_axcat = cp_axcat[:, ~cp_axcat.var.gene_ids.duplicated()]
    cp_axcat.obs_names_make_unique()
    cp_axcat.write(pathlib.PurePath(cp_axcat_norm_file))
    #
    toc = time.perf_counter()
    print(f"PREP (inc. CONCAT) in {toc - tic:0.4f} seconds")
    tic = time.perf_counter()
    #
    # axcat = ad.concat(au_objects)
    scement.combat(cp_axcat, key='batch', covariates=['age', 'stress'], inplace=True)
    sc.tl.pca(cp_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(cp_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(cp_axcat)
    sc.tl.umap(cp_axcat)
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_plot_file)
    cp_axcat.write(pathlib.PurePath(cp_axcat_results_file))
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")


def main(rtype):
    if rtype == "2ds":
        run_wt_gt_2dsets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "4ds":
        run_wt_gt_4dsets(DATA_DIR, ANALYSIS_OUT_DIR)
    elif rtype == "col0":
        run_col0_dsets(DATA_DIR, ANALYSIS_OUT_DIR)
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
    rtypes = ["2ds", "4ds", "col0", "copilot", "jb_gala", "gala", "default"]
    if len(sys.argv) <= 1:
        print("Usage: " + sys.argv[0] + " run_name [" + "/".join(rtypes) + "]")
    else:
        rtype = sys.argv[1]
        main(rtype)
