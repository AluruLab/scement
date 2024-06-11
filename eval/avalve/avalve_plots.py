import scanpy as sc
import anndata as an
import pathlib

# Assuming the analyses results are available in ANALYSIS_RESULTS_DIR
# Generate UMAP plots in PLOTS_OUT_DIR
ANALYSIS_RESULTS_DIR = "./avalve/"
PLOTS_OUT_DIR = "./avalve/"
COMBAT_RESULTS_FILE = "combat_avalve"
COMBAT_UMAP_FILE = "combat_avalve_renamed.png"
SCANORMA_RESULTS_FILE = "scanorama_avalve"
SCANORMA_UMAP_FILE = "scanorama_avalve_renamed.png"


def main():
    #
    sc._settings.ScanpyConfig.figdir = pathlib.Path(PLOTS_OUT_DIR)
    h5ad_file = ANALYSIS_RESULTS_DIR + COMBAT_RESULTS_FILE + ".h5ad"
    axcat = an.read_h5ad(h5ad_file)
    axcat = axcat[axcat.obs['cell_type'] != "Unknown"]
    axcat = axcat[axcat.obs["cell_type"] != "Monocytes"]
    axcat.obs["batch"] = axcat.obs["batch"].cat.rename_categories(  # type:ignore
        {"H1": "Healthy-1", "H2": "Healthy-2",
         "C3": "Diseased-1", "C4": "Diseased-2"})
    axcat.obs["cell_stype"] = axcat.obs["cell_stype"].cat.rename_categories(  # type:ignore
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
    sc.pl.umap(axcat, color=['batch', 'cell_stype'], save=COMBAT_UMAP_FILE)
    #
    sc._settings.ScanpyConfig.figdir = pathlib.Path(PLOTS_OUT_DIR)
    axcat = an.read_h5ad(ANALYSIS_RESULTS_DIR + SCANORMA_RESULTS_FILE + ".h5ad")
    axcat = axcat[axcat.obs['cell_type'] != "Unknown"]
    axcat = axcat[axcat.obs["cell_type"] != "Monocytes"]
    axcat.obs["batch"] = axcat.obs["batch"].cat.rename_categories(  # type:ignore
        {"H1": "Healthy-1", "H2": "Healthy-2",
         "C3": "Diseased-1", "C4": "Diseased-2"})
    axcat.obs["cell_stype"] = axcat.obs["cell_stype"].cat.rename_categories(  # type:ignore
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
    sc.pl.umap(axcat, color=['batch', 'cell_stype'], save=SCANORMA_UMAP_FILE)


if __name__ == "__main__":
    main()
