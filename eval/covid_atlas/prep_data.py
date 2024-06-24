import anndata as an
import scanpy as sc
import pandas as pd
import sys
import time
import argparse


def load_mtx(mat_dir):
    adata = sc.read_mtx(mat_dir + '/matrix.mtx.gz')
    adata_bc = pd.read_csv(mat_dir + '/barcodes.tsv.gz', sep='\t', header=None,
                           names=['cell_id'])  # type: ignore
    adata_features = pd.read_csv(mat_dir + '/features.tsv.gz',
                                 header=None, sep='\t',  # type: ignore
                                 names=['gene_id', 'gene_name'])  # type: ignore
    adata = adata.T
    adata.obs = adata_bc
    adata.obs.index = adata.obs['cell_id']  # type: ignore
    adata.var = adata_features
    adata.var.index = adata.var['gene_id']  # type: ignore
    return adata


def main(mat_dir, output_dir):
    tic = time.perf_counter()
    adx = load_mtx(mat_dir)
    print(adx, adx.var.shape, adx.obs.shape)
    toc = time.perf_counter()
    print(f"Completed matrix loading in {toc - tic:0.4f} seconds")
    #
    mtdx = pd.read_csv("./cell_annotation.csv")
    tic = time.perf_counter()
    for x in list(set(mtdx['sampleID'])):
        cnames = mtdx.loc[mtdx['sampleID'] == x, 'cellName']
        out_file = output_dir + "/" + x + ".h5ad"
        sample_adx = adx[cnames, ]
        sample_adx.write_h5ad(out_file)
    toc = time.perf_counter()
    print(f"Write matrices in {toc - tic:0.4f} seconds")
    #
    tic = time.perf_counter()
    for x in list(set(mtdx['sampleID'])):
        hd_file = output_dir + "/" + x + ".h5ad"
        sample_adx = an.read_h5ad(hd_file)
        print(sample_adx.X.shape, sample_adx.X.nnz)  # type: ignore
    toc = time.perf_counter()
    print(f"Check matrices in {toc - tic:0.4f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split COVID mtx as AnnData objects')
    help_str = "Path to the directory with matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz"
    parser.add_argument("matrix_dir", type=str, help=help_str)
    help_str = "Path to the output directory"
    parser.add_argument("output_dir", type=str, help=help_str)
    in_args = parser.parse_args(sys.argv[1:])
    print(in_args)
    print("Dir : ", in_args.matrix_dir, in_args.output_dir)
    main(in_args.matrix_dir, in_args.output_dir)
