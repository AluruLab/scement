import anndata as an
import pandas as pd
import scib

ATH_OUTPUT_DIR = "./athaliana/athaliana2"
ATH2_OUTPUT_DIR = "./athaliana/athaliana"
AVALVE_OUTPUT_DIR = "./avalve/avalve/"


def load_seurat_ad(seurat_ad, combat_ad):
    seurat_ad.obs["cell_type"] = combat_ad.obs.loc[seurat_ad.obs.index, "cell_type"]
    seurat_ad.obs["batch"] = combat_ad.obs.loc[seurat_ad.obs.index, "batch"]
    seurat_ad.var["gene_ids"] = seurat_ad.var["features"]
    if "cell_stype" in combat_ad.obs.columns:
        seurat_ad.obs["cell_stype"] = combat_ad.obs.loc[seurat_ad.obs.index,
                                                        "cell_stype"]
    if "cell_subtype" in combat_ad.obs.columns:
        seurat_ad.obs["cell_subtype"] = combat_ad.obs.loc[seurat_ad.obs.index,
                                                          "cell_subtype"]
    return seurat_ad


def load_fastmnn_ad(seurat_ad, combat_ad):
    # seurat_ad = an.read_h5ad(seurat_h5ad_file)
    seurat_ad.obs["cell_type"] = combat_ad.obs.loc[seurat_ad.obs.index, "cell_type"]
    if "cell_stype" in combat_ad.obs.columns:
        seurat_ad.obs["cell_stype"] = combat_ad.obs.loc[
            seurat_ad.obs.index, "cell_stype"
        ]
    if "cell_subtype" in combat_ad.obs.columns:
        seurat_ad.obs["cell_subtype"] = combat_ad.obs.loc[
            seurat_ad.obs.index, "cell_subtype"
        ]
    seurat_ad.obs["batch"] = combat_ad.obs.loc[seurat_ad.obs.index, "batch"]
    seurat_ad.var["gene_ids"] = seurat_ad.var["features"]
    seurat_ad.obsm["X_pca"] = seurat_ad.obsm["X_corrected"]
    return seurat_ad


def fastmnn_ath_metrics(output_dir="./athaliana"):
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala.h5ad"
    cp_axcat_results_file = output_dir + "/" + "fastmnn_jb_gala.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_fastmnn_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx4 = pd.DataFrame(mx)
    print(dmx4)
    rcx = {dmx4.columns[0]: "SEURAT-ATH-JB-GALA"}  # type:ignore
    dmx4.rename(columns=rcx, inplace=True)  # type:ignore
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    cp_axcat_results_file = output_dir + "/" + "fastmnn_jb_gala_flt2.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_fastmnn_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    print(dmx5)
    rcx = {dmx5.columns[0]: "SEURAT-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = pd.concat([dmx4, dmx5], axis=1)
    print(fdx)
    fdx.to_csv("athaliana/fastmnn_metrics.tsv", sep="\t")
    return dmx5


def seurat_ath_metrics(output_dir="./athaliana"):
    #
    # JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala.h5ad"
    cp_axcat_results_file = output_dir + "/seurat_ath_jb_gala.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx2 = pd.DataFrame(mx)
    rcx = {dmx2.columns[0]: "SEURAT-ATH-JB-GALA"}  # type:ignore
    dmx2.rename(columns=rcx, inplace=True)  # type:ignore
    print(dmx2)
    #  Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_gala.h5ad"
    cp_axcat_results_file = output_dir + "/" + "seurat_ath_gala.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx3 = pd.DataFrame(mx)
    rcx = {dmx3.columns[0]: "SEURAT-ATH-GALA"}  # type:ignore
    dmx3.rename(columns=rcx, inplace=True)  # type:ignore
    # Filter1 JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt1.h5ad"
    cp_axcat_results_file = output_dir + "/" + "seurat_ath_jb_gala_flt1.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx4 = pd.DataFrame(mx)
    rcx = {dmx4.columns[0]: "SEURAT-ATH-FLT1-JB-GALA"}  # type:ignore
    dmx4.rename(columns=rcx, inplace=True)  # type:ignore
    # Filter2 JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    cp_axcat_results_file = output_dir + "/" + "seurat_ath_jb_gala_flt2.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    rcx = {dmx5.columns[0]: "SEURAT-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = pd.concat([dmx2, dmx3, dmx4, dmx5], axis=1)
    print(fdx)
    fdx.to_csv("athaliana/seurat_metrics.tsv", sep="\t")
    return fdx


def scanorama_ath_metrics(output_dir="./athaliana/"):
    # Copilot
    #  cp_axcat_norm_file = output_dir + "/" + "scanpy_ath_coplit_norm.h5ad"
    #  file_sfx = "_ath_copilot"
    #  cp_axcat_results_file = output_dir + "/scanorama" + file_sfx + ".h5ad"
    #  adata = an.read_h5ad(cp_axcat_norm_file)
    #  adata_int = an.read_h5ad(cp_axcat_results_file)
    #  mx = scib.metrics.metrics(adata, adata_int, 'batch', 'cell_type',
    #                              ari_=True, nmi_=True, isolated_labels_asw_=True,
    #                              hvg_score_=True, pcr_=True, silhouette_=True,
    #                              isolated_labels_=True, graph_conn_=True)
    #  dmx1 = pd.DataFrame(mx)
    #  rcx = {dmx1.columns[0]: "SCANORAMA-ATH-COPILOT"} #type:ignore
    #  dmx1.rename(columns=rcx, inplace=True) #type:ignore
    # JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala.h5ad"
    file_sfx = "ath_jb_gala"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        embed="X_scanorama",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx2 = pd.DataFrame(mx)
    rcx = {dmx2.columns[0]: "SCANORAMA-ATH-JB-GALA"}  # type:ignore
    dmx2.rename(columns=rcx, inplace=True)  # type:ignore
    # print(dmx2)
    # Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_gala.h5ad"
    file_sfx = "ath_gala"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    # mx = scib.metrics.metrics(adata, adata_int, 'batch', 'cell_type',
    #                          embed='X_scanorama',
    #                          ari_=True, nmi_=True, isolated_labels_asw_=True,
    #                          hvg_score_=True, pcr_=True, silhouette_=True,
    #                          isolated_labels_=True, graph_conn_=True)
    # dmx3 = pd.DataFrame(mx)
    # rcx = {dmx3.columns[0]: "SCANORAMA-ATH-GALA"} #type:ignore
    # dmx3.rename(columns=rcx, inplace=True) #type:ignore
    dmx3 = dmx2
    # Filter1 JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt1.h5ad"
    file_sfx = "ath_jb_gala_flt1"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        embed="X_scanorama",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx4 = pd.DataFrame(mx)
    rcx = {dmx4.columns[0]: "SCANORAMA-ATH-FLT1-JB-GALA"}  # type:ignore
    dmx4.rename(columns=rcx, inplace=True)  # type:ignore
    # Filter2 JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    file_sfx = "ath_jb_gala_flt2"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        embed="X_scanorama",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    rcx = {dmx5.columns[0]: "SCANORAMA-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = pd.concat([dmx2, dmx3, dmx4, dmx5], axis=1)
    print(fdx)
    fdx.to_csv("athaliana/scanorama_metrics.tsv", sep="\t")
    return fdx


def combat_ath_metrics(output_dir="./athaliana/"):
    # copilot
    #  cp_axcat_results_file = output_dir + "/" + "scanpy_ath_copilot.h5ad"
    #  cp_axcat_norm_file = output_dir + "/" + "scanpy_ath_coplit_norm.h5ad"
    #  adata = an.read_h5ad(cp_axcat_norm_file)
    #  adata_int = an.read_h5ad(cp_axcat_results_file)
    #  mx = scib.metrics.metrics(adata, adata_int, 'batch', 'cell_type',
    #                             ari_=True, nmi_=True, isolated_labels_asw_=True,
    #                             hvg_score_=True, pcr_=True, silhouette_=True,
    #                             isolated_labels_=True, graph_conn_=True)
    #  dmx1 = pd.DataFrame(mx)
    #  rcx = {dmx1.columns[0]: "COMBAT-ATH-COPILOT"} #type:ignore
    #  dmx1.rename(columns=rcx, inplace=True) #type:ignore
    # JB Gala
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx2 = pd.DataFrame(mx)
    rcx = {dmx2.columns[0]: "COMBAT-ATH-JB-GALA"}  # type:ignore
    dmx2.rename(columns=rcx, inplace=True)  # type:ignore
    # Gala
    cp_axcat_results_file = output_dir + "/" + "combat_ath_gala.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_gala.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx3 = pd.DataFrame(mx)
    rcx = {dmx3.columns[0]: "COMBAT-ATH-GALA"}  # type:ignore
    dmx3.rename(columns=rcx, inplace=True)  # type:ignore
    # Filter1 JB Gala
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala_flt1.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt1.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx4 = pd.DataFrame(mx)
    rcx = {dmx4.columns[0]: "COMBAT-ATH-FLT1-JB-GALA"}  # type:ignore
    dmx4.rename(columns=rcx, inplace=True)  # type:ignore
    # Filter2 JB Gala
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala_flt2.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    rcx = {dmx5.columns[0]: "COMBAT-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = pd.concat([dmx2, dmx3, dmx4, dmx5], axis=1)
    print(fdx)
    fdx.to_csv("athaliana/combat_metrics.tsv", sep="\t")
    return fdx


def fastmnn_ath_metrics2(output_dir="./athaliana2"):
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    cp_axcat_results_file = output_dir + "/" + "fastmnn_jb_gala_flt2.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_fastmnn_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    print(dmx5)
    rcx = {dmx5.columns[0]: "SEURAT-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = dmx5
    print(fdx)
    fdx.to_csv("athaliana2/fastmnn_metrics.tsv", sep="\t")
    return dmx5


def seurat_ath_metrics2(output_dir="./athaliana"):
    # Filter2 JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    cp_axcat_results_file = output_dir + "/" + "seurat_ath_jb_gala_flt2.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    rcx = {dmx5.columns[0]: "SEURAT-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = dmx5
    print(fdx)
    fdx.to_csv("athaliana2/seurat_metrics.tsv", sep="\t")
    return fdx


def scanorama_ath_metrics2(output_dir="./athaliana/"):
    # Filter2 JB Gala
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    file_sfx = "ath_jb_gala_flt2"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        embed="X_scanorama",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    rcx = {dmx5.columns[0]: "SCANORAMA-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = dmx5
    print(fdx)
    fdx.to_csv("athaliana2/scanorama_metrics.tsv", sep="\t")
    return fdx


def combat_ath_metrics2(output_dir="./athaliana/"):
    # Filter2 JB Gala
    cp_axcat_results_file = output_dir + "/" + "combat_ath_jb_gala_flt2.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_ath_jb_gala_flt2.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    rcx = {dmx5.columns[0]: "COMBAT-ATH-FLT2-JB-GALA"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    fdx = dmx5
    print(fdx)
    fdx.to_csv("athaliana2/combat_metrics.tsv", sep="\t")
    return fdx


def ath_metrics():
    fx4 = fastmnn_ath_metrics(output_dir=ATH_OUTPUT_DIR)
    fx1 = combat_ath_metrics(output_dir=ATH_OUTPUT_DIR)
    fx2 = scanorama_ath_metrics(output_dir=ATH_OUTPUT_DIR)
    fx3 = seurat_ath_metrics(output_dir=ATH_OUTPUT_DIR)
    fdx = pd.concat([fx1, fx2, fx3, fx4], axis=1)
    fdx.to_csv("athaliana/scib_metrics.tsv", sep="\t")


def ath_metrics2():
    fx4 = fastmnn_ath_metrics2(output_dir=ATH2_OUTPUT_DIR)
    fx1 = combat_ath_metrics2(output_dir=ATH_OUTPUT_DIR)
    fx2 = scanorama_ath_metrics2(output_dir=ATH_OUTPUT_DIR)
    fx3 = seurat_ath_metrics2(output_dir=ATH_OUTPUT_DIR)
    fdx = pd.concat([fx1, fx2, fx3, fx4], axis=1)
    fdx.to_csv(ATH2_OUTPUT_DIR + "/scib_metrics.tsv", sep="\t")


def avalve_metrics(output_dir):
    cp_axcat_norm_file = output_dir + "2/" + "scanpy_normalized_avalve.h5ad"
    cp_axcat_results_file = output_dir + "2/seurat_avalve.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_type",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx2 = pd.DataFrame(mx)
    rcx = {dmx2.columns[0]: "SEURAT-AVALVE"}  # type:ignore
    dmx2.rename(columns=rcx, inplace=True)  # type:ignore
    print(dmx2)
    dmx2.to_csv(output_dir + "/seurat_metrics2.tsv", sep="\t")
    #
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_avalve.h5ad"
    file_sfx = "avalve"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_stype",
        embed="X_scanorama",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx3 = pd.DataFrame(mx)
    rcx = {dmx3.columns[0]: "SCANORAMA-AVALVE"}  # type:ignore
    dmx3.rename(columns=rcx, inplace=True)  # type:ignore
    print(dmx3)
    dmx3.to_csv(output_dir + "/scanorama_metrics2.tsv", sep="\t")
    #
    cp_axcat_results_file = output_dir + "/" + "combat_avalve.h5ad"
    cp_axcat_norm_file = output_dir + "/" + "scanpy_normalized_avalve.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_stype",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx4 = pd.DataFrame(mx)
    rcx = {dmx4.columns[0]: "COMBAT-AVALVE"}  # type:ignore
    dmx4.rename(columns=rcx, inplace=True)  # type:ignore
    dmx4.to_csv(AVALVE_OUTPUT_DIR + "/combat_metrics2.tsv", sep="\t")
    print(dmx4)
    cp_axcat_norm_file = output_dir + "/" + "combat_avalve.h5ad"
    cp_axcat_results_file = output_dir + "/" + "fastmnn_avalve.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata_int = load_fastmnn_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_stype",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=True,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    print(dmx5)
    rcx = {dmx5.columns[0]: "FASTMNN-AVALVE"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    dmx5.to_csv(AVALVE_OUTPUT_DIR + "/fastmnn_metrics2.tsv", sep="\t")
    return dmx5


def pbmc_metrics(output_dir="./pbmc"):
    cp_axcat_norm_file = output_dir + "/" + "scanpy_merge_pbmc.h5ad"
    cp_axcat_results_file = output_dir + "/seurat_pbmc.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata = adata[adata.obs.index.intersection(adata_int.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index.intersection(adata.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index,  # type:ignore
                          adata_int.var.index.intersection(adata.var.index)]
    adata = adata[adata.obs.index, adata.var.index.intersection(
        adata_int.var.index)]  # type:ignore
    adata_int = load_seurat_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_subtype",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=False,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx2 = pd.DataFrame(mx)
    rcx = {dmx2.columns[0]: "SEURAT-PBMC"}  # type:ignore
    dmx2.rename(columns=rcx, inplace=True)  # type:ignore
    print(dmx2)
    dmx2.to_csv("pbmc/seurat_metrics2.tsv", sep="\t")
    #
    cp_axcat_norm_file = output_dir + "/" + "scanpy_merge_pbmc.h5ad"
    file_sfx = "pbmc"
    cp_axcat_results_file = output_dir + "/scanorama_" + file_sfx + ".h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata = adata[adata.obs.index.intersection(adata_int.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index.intersection(adata.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index,
                          adata_int.var.index.intersection(adata.var.index)]  # type:ignore
    adata = adata[adata.obs.index,
                  adata.var.index.intersection(adata_int.var.index)]  # type:ignore
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_subtype",
        embed="X_scanorama",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=False,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx3 = pd.DataFrame(mx)
    rcx = {dmx3.columns[0]: "SCANORAMA-PBMC"}  # type:ignore
    dmx3.rename(columns=rcx, inplace=True)  # type:ignore
    print(dmx3)
    dmx3.to_csv("pbmc/scanorama_metrics2.tsv", sep="\t")
    #
    cp_axcat_norm_file = output_dir + "/" + "scanpy_merge_pbmc.h5ad"
    cp_axcat_results_file = output_dir + "/" + "combat_pbmc.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata = adata[adata.obs.index.intersection(adata_int.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index.intersection(adata.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index,  # type:ignore
                          adata_int.var.index.intersection(adata.var.index)]
    adata = adata[adata.obs.index,  # type:ignore
                  adata.var.index.intersection(adata_int.var.index)]
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_subtype",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=False,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx4 = pd.DataFrame(mx)
    rcx = {dmx4.columns[0]: "COMBAT-PBMC"}  # type:ignore
    dmx4.rename(columns=rcx, inplace=True)  # type:ignore
    dmx4.to_csv("pbmc/combat_metrics2.tsv", sep="\t")
    print(dmx4)
    cp_axcat_norm_file = output_dir + "/" + "combat_pbmc.h5ad"
    cp_axcat_results_file = output_dir + "/" + "fastmnn_pbmc.h5ad"
    adata = an.read_h5ad(cp_axcat_norm_file)
    adata_int = an.read_h5ad(cp_axcat_results_file)
    adata = adata[adata.obs.index.intersection(adata_int.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index.intersection(adata.obs.index), ]  # type:ignore
    adata_int = adata_int[adata_int.obs.index,  # type:ignore
                          adata_int.var.index.intersection(adata.var.index)]
    adata = adata[adata.obs.index,  # type:ignore
                  adata.var.index.intersection(adata_int.var.index)]
    adata_int = load_fastmnn_ad(adata_int, adata)
    mx = scib.metrics.metrics(
        adata,
        adata_int,
        "batch",
        "cell_subtype",
        ari_=True,
        nmi_=True,
        isolated_labels_asw_=True,
        hvg_score_=True,
        pcr_=True,
        silhouette_=False,
        isolated_labels_=True,
        graph_conn_=True,
    )
    dmx5 = pd.DataFrame(mx)
    print(dmx5)
    rcx = {dmx5.columns[0]: "FASTMNN-PBMC"}  # type:ignore
    dmx5.rename(columns=rcx, inplace=True)  # type:ignore
    dmx5.to_csv("pbmc/fastmnn_metrics2.tsv", sep="\t")

    # fdx = pd.concat([dmx2, dmx3, dmx4, dmx5], axis=1)
    # print(fdx)
    # fdx.to_csv("pbmc/scib_metrics.tsv", sep="\t")
    return dmx5


def main():
    # ath_metrics(ATH_OUTPUT_DIR)
    ath_metrics2()
    avalve_metrics(AVALVE_OUTPUT_DIR)
    # pbmc_metrics()


if __name__ == "__main__":
    main()
