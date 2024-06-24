import scanpy as sc
import anndata as an
import sys
import pathlib

ATH_DATA_DIR = "./athaliana/athaliana/"
AVALVE_DATA_DIR = "./avalve/avalve/"


def generate_plot(seurat_h5ad_file, combat_h5ad_file, output_dir, file_sfx,
                  filter_unknown=False):
    seurat_ad = an.read_h5ad(seurat_h5ad_file)
    combat_ad = an.read_h5ad(combat_h5ad_file)
    seurat_ad = seurat_ad[seurat_ad.obs.index.intersection(combat_ad.obs.index), ]  # type:ignore
    seurat_ad.obs["cell_type"] = combat_ad.obs.loc[seurat_ad.obs.index,  # type:ignore
                                                   "cell_type"]
    seurat_ad.obs["batch"] = combat_ad.obs.loc[seurat_ad.obs.index, "batch"]
    seurat_ad.var["gene_ids"] = seurat_ad.var["features"]
    if "cell_stype" in combat_ad.obs.columns:
        seurat_ad.obs["cell_stype"] = combat_ad.obs.loc[seurat_ad.obs.index,
                                                        "cell_stype"]
    # seurat_ad.obs["Assay"] = combat_ad.obs.loc[seurat_ad.obs.index, "Assay"]
    if "health" in combat_ad.obs.columns:
        seurat_ad.obs["health"] = combat_ad.obs.loc[seurat_ad.obs.index, "health"]
        seurat_ad.obs["condition"] = seurat_ad.obs["health"].cat.rename_categories(  # type:ignore
            {"H": "Health", "C": "Calcified"})
        seurat_ad.obs["batch"] = seurat_ad.obs["batch"].cat.rename_categories({  # type:ignore
            "H1": "Healthy-1",
            "H2": "Healthy-2",
            "C3": "Diseased-1",
            "C4": "Diseased-2",
        })
        seurat_ad.obs["cell_stype"] = seurat_ad.obs["cell_stype"].cat.rename_categories({  # type:ignore
            "CLymphocytes": "Diseased-Lymphocytes",
            "CMacrophages": "Diseased-Macrophages",
            "CVEC": "Diseased-VEC",
            "CVIC": "Diseased-VIC",
            "HLymphocytes": "Healthy-Lymphocytes",
            "HMacrophages": "Healthy-Macrophages",
            "HVEC": "Healthy-VEC",
            "HVIC": "Healthy-VIC",
        })
    if "cell_subtype" in combat_ad.obs.columns:
        seurat_ad.obs["cell_stype"] = combat_ad.obs.loc[seurat_ad.obs.index,
                                                        "cell_subtype"]
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    if filter_unknown:
        plots = ["batch", "cell_type"]
        plots_sfilter = ["batch", "cell_stype"]
        filter_seurat_plot_file = "seurat_" + file_sfx + ".png"
        seurat_ad = seurat_ad[seurat_ad.obs["cell_type"] != "Unknown"]
        sc.pl.umap(seurat_ad, color=plots, save=filter_seurat_plot_file)
        filter_seurat_plot_file = "seurat_" + file_sfx + "_sfilter.png"
        sc.pl.umap(seurat_ad, color=plots_sfilter, save=filter_seurat_plot_file)
        seurat_ad2 = seurat_ad[seurat_ad.obs["cell_type"] != "Monocytes"]
        filter2_seurat_plot_file = "seurat_" + file_sfx + "_final.png"
        sc.pl.umap(seurat_ad2, color=plots_sfilter,
                   save=filter2_seurat_plot_file)
    else:
        if "cell_stype" in seurat_ad.obs.columns:
            cp_types = ["batch", "cell_type", "cell_stype"]
        else:
            cp_types = ["batch", "cell_type"]
        print(cp_types)
        seurat_plot_file = "seurat_" + file_sfx + ".png"
        sc.pl.umap(seurat_ad, color=cp_types, save=seurat_plot_file)
    return seurat_ad


def main(rtype):
    if rtype == "jb_gala":
        seurat_h5ad_file = ATH_DATA_DIR + "/seurat_ath_jb_gala.h5ad"
        combat_h5ad_file = ATH_DATA_DIR + "/scanpy_normalized_ath_jb_gala.h5ad"
        return generate_plot(seurat_h5ad_file, combat_h5ad_file, ATH_DATA_DIR,
                             "ath_jb_gala", False)
    elif rtype == "gala":
        seurat_h5ad_file = ATH_DATA_DIR + "/seurat_ath_gala.h5ad"
        combat_h5ad_file = ATH_DATA_DIR + "/scanpy_normalized_ath_gala.h5ad"
        return generate_plot(seurat_h5ad_file, combat_h5ad_file, ATH_DATA_DIR,
                             "ath_gala", False)
    elif rtype == "copilot":
        seurat_h5ad_file = ATH_DATA_DIR + "/seurat_ath_copilot.h5ad"
        combat_h5ad_file = ATH_DATA_DIR + "/scanpy_normalized_ath_copilot.h5ad"
        return generate_plot(seurat_h5ad_file, combat_h5ad_file, ATH_DATA_DIR,
                             "ath_copilot", False)
    elif rtype == "jb_gala_flt1":
        seurat_h5ad_file = ATH_DATA_DIR + "/seurat_ath_jb_gala_flt1.h5ad"
        combat_h5ad_file = ATH_DATA_DIR + "/scanpy_normalized_ath_jb_gala_flt1.h5ad"
        return generate_plot(seurat_h5ad_file, combat_h5ad_file, ATH_DATA_DIR,
                             "ath_jb_gala_flt1", False)
    elif rtype == "jb_gala_flt2":
        seurat_h5ad_file = ATH_DATA_DIR + "/seurat_ath_jb_gala_flt2.h5ad"
        combat_h5ad_file = ATH_DATA_DIR + "/scanpy_normalized_ath_jb_gala_flt2.h5ad"
        return generate_plot(seurat_h5ad_file, combat_h5ad_file, ATH_DATA_DIR,
                             "ath_jb_gala_flt2", False)
    elif rtype == "avalve":
        seurat_h5ad_file = AVALVE_DATA_DIR + "/seurat_avalve.h5ad"
        combat_h5ad_file = AVALVE_DATA_DIR + "/scanpy_normalized_avalve.h5ad"
        return generate_plot(seurat_h5ad_file, combat_h5ad_file,
                             AVALVE_DATA_DIR, "avalve", True)


if __name__ == "__main__":
    rtypes = ["copilot", "jb_gala", "gala", "jb_gala_flt1", "avalve"]
    if len(sys.argv) <= 1:
        print("Usage: " + sys.argv[0] + " run_name [" + "/".join(rtypes) + "]")
    else:
        rtype = sys.argv[1]
        main(rtype)
