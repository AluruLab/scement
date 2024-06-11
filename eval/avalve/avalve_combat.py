import scanpy as sc
import anndata as an
import time
import pathlib

# Generate analyses with COMBAT
# source data is from DATA_DIR, and output is ANALYSIS_RESULTS_DIR
DATA_DIR = "./aortic_valve/"
ANALYSIS_RESULTS_DIR = "./avalve/"


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
    #
    # axcat = ad.concat(au_objects)
    axcat = an.concat(au_objects, merge='same')
    axcat.obs_names_make_unique()
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    return axcat


def run_avalve_dsets1(data_dir, output_dir):
    batches = ['H1', 'H2', 'C3', 'C4']
    #
    tic = time.perf_counter()
    full_ads = [an.read_h5ad(data_dir + bx + ".h5ad") for bx in batches]
    for ix in range(len(full_ads)):
        bx = batches[ix]
        full_ads[ix].obs['health'] = bx[0]
        full_ads[ix].obs['cell_stype'] = full_ads[ix].obs['health'].astype(  # type:ignore
            'str') + full_ads[ix].obs['cell_type'].astype('str')  # type:ignore
        filter_scale(full_ads[ix])
    cp_axcat = prep_combat_intersect(full_ads)
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    cp_axcat_plot_file = "combat_avalve.png"
    cp_axcat_results_file = output_dir + "/" + "combat_avalve.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_avalve.h5ad"
    cp_axcat = cp_axcat[:, ~cp_axcat.var.gene_ids.duplicated()]
    cp_axcat.obs_names_make_unique()
    cp_axcat.write(cp_axcat_norm_file)
    toc = time.perf_counter()
    print(f"TOTAL PREP (w. INTERSECT) in {toc - tic:0.4f} seconds")
    #
    tic = time.perf_counter()
    sc.pp.combat(cp_axcat, key='batch', inplace=True)
    sc.tl.pca(cp_axcat, svd_solver='arpack', n_comps=20)
    sc.pp.neighbors(cp_axcat, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(cp_axcat)
    sc.tl.umap(cp_axcat)
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_plot_file)
    cp_axcat.write(cp_axcat_results_file)
    toc = time.perf_counter()
    print(f"COMBAT in {toc - tic:0.4f} seconds")
    #
    cp_axcat_filter_plot_file = "combat_avalve_filter.png"
    cp_axcat = cp_axcat[cp_axcat.obs['cell_type'] != "Unknown"]
    sc.pl.umap(cp_axcat, color=['batch', 'cell_type'], save=cp_axcat_filter_plot_file)
    cp_axcat_sfilter_plot_file = "combat_avalve_sfilter.png"
    sc.pl.umap(cp_axcat, color=['batch', 'cell_stype'], save=cp_axcat_sfilter_plot_file)
    #
    cp_axcat2 = cp_axcat[cp_axcat.obs["cell_type"] != "Monocytes"]
    cp_axcat2.obs["batch"] = cp_axcat2.obs["batch"].cat.rename_categories(  # type:ignore
        {"H1": "Healthy-1", "H2": "Healthy-2",
         "C3": "Diseased-1", "C4": "Diseased-2"})
    cp_axcat2.obs["cell_stype"] = cp_axcat2.obs["cell_stype"].cat.rename_categories(  # type:ignore
        {
            "CLymphocytes": "Diseased-Lymphocytes",
            "CMacrophages": "Diseased-Macrophages",
            "CVEC": "Diseased-VEC",
            "CVIC": "Diseased-VIC",
            "HLymphocytes": "Healthy-Lymphocytes",
            "HMacrophages": "Healthy-Macrophages",
            "HVEC": "Healthy-VEC",
            "HVIC": "Healthy-VIC"
        })
    cp_axcat_sfilter2_plot_file = "combat_avalve_final.png"
    sc.pl.umap(cp_axcat2, color=['batch', 'cell_stype'], save=cp_axcat_sfilter2_plot_file)


def main():
    run_avalve_dsets1(DATA_DIR, ANALYSIS_RESULTS_DIR)


if __name__ == "__main__":
    main()
