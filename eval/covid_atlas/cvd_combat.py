import pandas as pd
import scanpy as sc
import anndata as an
import sys
import time
import scement as sct
import pathlib
import argparse


def load_metadata(nmin=1000):
    dfx = pd.read_excel("meta/covid_metadata.xlsx",  # type: ignore
                        index_col="sampleID")   # type: ignore
    filter_types = ["fresh BALF", "fresh Sputum", "fresh PFMC"]
    dfx = dfx.loc[~dfx["Sample type"].isin(filter_types), ]
    # vcx = dfx["Sample type"].value_counts()
    #
    qcx = pd.read_excel("meta/covid_metadata.xlsx", sheet_name="qc",
                        index_col="sampleID")  # type:ignore
    qcx = qcx.loc[dfx.index]
    qcx = qcx.loc[qcx["Ncells_manual_removal"] > nmin]
    dfx = dfx.loc[qcx.index]
    cellfx = pd.read_csv("meta/cell_annotation.csv.gz", index_col='cellName')
    return qcx, dfx, cellfx


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


def load_filter_data(hd_dir, ctrl_kp, dis_kp, ngenes):
    _, dfx, cellfx = load_metadata()
    tic = time.perf_counter()
    ctrl_ads = [an.read_h5ad(hd_dir + "/" + x + ".h5ad") for x in ctrl_kp]
    for i, x in enumerate(ctrl_kp):
        ctrl_ads[i].obs["pid"] = x
        ctrl_ads[i].obs["batch"] = dfx.loc[x]["Batches"]
        ctrl_ads[i].obs["covid"] = dfx.loc[x]["SARS-CoV-2"]
        ctrl_ads[i].obs["severity"] = dfx.loc[x]["CoVID-19 severity"]
        ctrl_ads[i].obs["type"] = dfx.loc[x]['Sample type']
        ctrl_ads[i].obs["cell_type"] = cellfx.loc[  # type:ignore
            ctrl_ads[i].obs.index, "majorType"]
        ctrl_ads[i].obs["cell_sub_type"] = cellfx.loc[  # type:ignore
            ctrl_ads[i].obs.index, "celltype"]
        filter_scale(ctrl_ads[i], ngenes)
    toc = time.perf_counter()
    print(f"CTRL LOAD in {toc - tic:0.4f} seconds")
    #
    tic = time.perf_counter()
    dis_ads = [an.read_h5ad(hd_dir + "/" + x + ".h5ad") for x in dis_kp]
    for i, x in enumerate(dis_kp):
        dis_ads[i].obs["pid"] = x
        dis_ads[i].obs["batch"] = dfx.loc[x]["Batches"]
        dis_ads[i].obs["covid"] = dfx.loc[x]["SARS-CoV-2"]
        dis_ads[i].obs["severity"] = dfx.loc[x]["CoVID-19 severity"]
        dis_ads[i].obs["type"] = dfx.loc[x]['Sample type']
        dis_ads[i].obs["cell_type"] = cellfx.loc[  # type:ignore
            dis_ads[i].obs.index, "majorType"]
        dis_ads[i].obs["cell_sub_type"] = cellfx.loc[  # type:ignore
            dis_ads[i].obs.index, "celltype"]
        filter_scale(dis_ads[i], ngenes)
    toc = time.perf_counter()
    print(f"DIS LOAD in {toc - tic:0.4f} seconds")
    #
    all_ads = ctrl_ads + dis_ads
    axcat = prep_combat_intersect(all_ads)
    return all_ads, axcat


def combat_and_plot(axcat, axcat_plot_file, axcat_file,
                    ntopg=5000, scement=True):
    if ntopg is not None:
        tic = time.perf_counter()
        sc.pp.highly_variable_genes(axcat, flavor="seurat", n_top_genes=ntopg)
        axcat = axcat[:, axcat.var.highly_variable]
        toc = time.perf_counter()
        print(f"HVG in {toc - tic:0.4f} seconds")
    print(axcat)
    tic = time.perf_counter()
    if scement:
        # tdc.combat(axcat, key='pid', covariates=['covid', 'type'],
        #            inplace=True)
        sct.sct_sparse(axcat, key='pid', covariates=['covid', 'type'],
                       inplace=True)
    else:
        sct.combat(axcat, key='pid', covariates=['covid', 'type'],
                   inplace=True)
        # sc.pp.combat(axcat, key='pid', covariates=['covid', 'type'],
        #              inplace=True)
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
    tic = time.perf_counter()
    sc.pl.umap(axcat, color=['pid', 'covid', 'cell_type'],
               save=axcat_plot_file)
    toc = time.perf_counter()
    print(f"PLOT in {toc - tic:0.4f} seconds")


def build_50k_ad(hd_dir, ngenes):
    ctrl_kp = ["S-HC025", "S-HC010", "S-HC013", "S-HC008", "S-HC019-2"]
    dis_kp = ["S-M035-1", "S-S035-2", "S-M061-2", "S-M028", "S-S079"]
    lsx , axcat = load_filter_data(hd_dir, ctrl_kp, dis_kp, ngenes)
    return lsx, axcat


def run_50K(hd_dir, gxtype, ngenes, output_dir=".",
            preprocess_only=False, scement=True):
    _, axcat = build_50k_ad(hd_dir, ngenes)
    print(axcat)
    pp_file = output_dir + "/" + "pp_50k_" + gxtype + ".h5ad"
    if scement:
        axcat_file = output_dir + "/" + "scement_50k_" + gxtype + ".h5ad"
        axcat_plot_file = "_scement_50k_" + gxtype + ".png"
    else:
        axcat_file = output_dir + "/" + "sc_combat_50k_" + gxtype + ".h5ad"
        axcat_plot_file = "_sc_combat_50k_" + gxtype + ".png"
    if preprocess_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file, ngenes, scement)


def build_320k_ad(hd_dir, ngenes):
    _, dfx, _ = load_metadata()
    ctrl_kp = list(dfx[dfx['SARS-CoV-2'] == 'negative'].index)
    dis_kp = ['S-M008-2', 'S-S062', 'S-M015', 'S-M004-3', 'S-S014',
              'S-M053', 'S-M061-2', 'S-M006', 'S-S083', 'S-S022-3',
              'S-M036-1', 'S-S016', 'S-S025', 'S-M046', 'S-S074-2',
              'S-S031', 'S-M010-2', 'S-S077', 'S-M009-5', 'S-M044-1',
              'S-S085-2', 'S-M008-1', 'S-M058-2', 'S-S016', 'S-M070',
              'S-M016', 'S-S030', 'S-S032-3']
    lsx, axcat = load_filter_data(hd_dir, ctrl_kp, dis_kp, ngenes)
    return lsx, axcat


def run_320K(hd_dir, gxtype, ngenes, output_dir=".",
             preprocess_only=False, scement=True):
    _, axcat = build_320k_ad(hd_dir, ngenes)
    print(axcat)
    pp_file = output_dir + "/" + "pp_320k_" + gxtype + ".h5ad"
    if scement:
        axcat_file = output_dir + "/" + "scement_320k_" + gxtype + ".h5ad"
        axcat_plot_file = "_scement_320k_" + gxtype + ".png"
    else:
        axcat_file = output_dir + "/" + "sc_combat_320k_" + gxtype + ".h5ad"
        axcat_plot_file = "_sc_combat_320k_" + gxtype + ".png"
    if preprocess_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file, ngenes, scement)


def build_500k_ad(hd_dir, ngenes):
    _, dfx, _ = load_metadata()
    ctrl_kp = list(dfx[dfx['SARS-CoV-2'] == 'negative'].index)
    dis_kp = ['S-S068', 'S-M018', 'S-S050', 'S-M066', 'S-S013', 'S-S035-1',
              'S-S001-2', 'S-S081', 'S-M037', 'S-M042-2', 'S-S078', 'S-S047',
              'S-M048', 'S-S018', 'S-M042-1', 'S-S050', 'S-M010-2', 'S-M009-2',
              'S-S084', 'S-S044', 'S-S047', 'S-M073', 'S-M050', 'S-M016',
              'S-M008-2', 'S-M078', 'S-M011', 'S-M054', 'S-S092', 'S-M007-4',
              'S-M064', 'S-S082', 'S-M004-1', 'S-S068', 'S-S035-4', 'S-M030',
              'S-S013', 'S-M004-6', 'S-S092', 'S-M033', 'S-S057', 'S-M041-1',
              'S-M063-1', 'S-M010-1', 'S-M069', 'S-M049', 'S-M059-2', 'S-S017',
              'S-S082', 'S-M005']
    lsx, axcat = load_filter_data(hd_dir, ctrl_kp, dis_kp, ngenes)
    return lsx, axcat


def run_500K(hd_dir, gxtype, ngenes, output_dir=".",
             preprocess_only=False, scement=True):
    _, axcat = build_500k_ad(hd_dir, ngenes)
    print(axcat)
    pp_file = output_dir + "/" + "pp_500k_" + gxtype + ".h5ad"
    if scement:
        axcat_file = output_dir + "/" + "scement_500k_" + gxtype + ".h5ad"
        axcat_plot_file = "_scement_500k_" + gxtype + ".png"
    else:
        axcat_file = output_dir + "/" + "sc_combat_500k_" + gxtype + ".h5ad"
        axcat_plot_file = "_sc_combat_500k_" + gxtype + ".png"
    if preprocess_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file, ngenes, scement)


def build_full_ad(hd_dir, ngenes):
    _, dfx, _ = load_metadata()
    ctrl_kp = list(dfx[dfx['SARS-CoV-2'] == 'negative'].index)
    dis_kp = list(dfx[dfx['SARS-CoV-2'] == 'positive'].index)
    lsx , axcat = load_filter_data(hd_dir, ctrl_kp, dis_kp, ngenes)
    return lsx, axcat


def run_full(hd_dir, gxtype, ngenes, output_dir=".",
             preprocess_only=False, scement=True):
    _, axcat = build_full_ad(hd_dir, ngenes)
    print(axcat)
    pp_file = output_dir + "/" + "pp_full_" + gxtype + ".h5ad"
    if scement:
        axcat_file = output_dir + "/" + "scement_full_" + gxtype + ".h5ad"
        axcat_plot_file = "_scement_full_" + gxtype + ".png"
    else:
        axcat_plot_file = "_sc_combat_full_" + gxtype + ".png"
        axcat_file = output_dir + "/" + "sc_combat_full_" + gxtype + ".h5ad"
    if preprocess_only:
        axcat.write_h5ad(pp_file)
    else:
        combat_and_plot(axcat, axcat_plot_file, axcat_file, ngenes, scement)


def main(hd_dir, cxtype, gxtype, ngenes, output_dir=".",
         preprocess_only=False, scement=True):
    sc._settings.ScanpyConfig.figdir = pathlib.Path(output_dir)
    if cxtype == "50K":
        run_50K(hd_dir, gxtype, ngenes, output_dir, preprocess_only, scement)
    elif cxtype == "320K":
        run_320K(hd_dir, gxtype, ngenes, output_dir, preprocess_only, scement)
    elif cxtype == "500K":
        run_500K(hd_dir, gxtype, ngenes, output_dir, preprocess_only, scement)
    elif cxtype == "full":
        run_full(hd_dir, gxtype, ngenes, output_dir, preprocess_only, scement)


if __name__ == "__main__":
    rctypes = ["50K", "320K", "500K", "full"]
    hd_dir = "./h5ad_full/"
    out_dir = "./"
    ngenes = {"500": 500, "1K": 1000, "5K": 5000, "10K": 10000, "full": None}
    rgtypes = list(ngenes.keys())
    parser = argparse.ArgumentParser(prog="CVD",
                                     description="CVD COMBAT",
                                     epilog="CVD HELP")
    parser.add_argument("cxtype", choices=rctypes)
    parser.add_argument("gxtype", choices=rgtypes)
    parser.add_argument("-p", "--preprocess_only", action="store_true")
    parser.add_argument("-s", "--scement", action="store_true")
    in_args = parser.parse_args(sys.argv[1:])
    print(hd_dir, in_args.cxtype, in_args.gxtype,
          ngenes[in_args.gxtype], out_dir, in_args.preprocess_only,
          in_args.scement)
    main(hd_dir, in_args.cxtype, in_args.gxtype,
         ngenes[in_args.gxtype], out_dir, in_args.preprocess_only, in_args.scement)
