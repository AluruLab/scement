import anndata as ad
import argparse
import json
import time


def main(json_data):
    data_dir = json_data['DATADIR']
    in_files = json_data['ADFILE']
    # batch_names = json_data['BATCH']
    out_dir = json_data['OUTDIR']
    out_file = json_data['OUTFILE']
    results_file = out_dir + "/" + out_file
    #
    tic = time.perf_counter()
    ad_objects = [ad.read_h5ad(data_dir + "/" + fx) for fx in in_files]
    for x in ad_objects:
        x.var.index = [y for y in x.var.gene_ids]
    au_objects = [ax[:, ~ax.var.gene_ids.duplicated()] for ax in ad_objects]
    toc = time.perf_counter()
    print(f"LOAD in {toc - tic:0.4f} seconds")
    print(ad_objects)
    print(au_objects)
    #
    tic = time.perf_counter()
    # axcat = ad.concat(au_objects)
    axcat = ad.concat(au_objects, merge='same')
    axcat.obs_names_make_unique()
    toc = time.perf_counter()
    print(f"CONCAT in {toc - tic:0.4f} seconds")
    print(axcat)
    #
    tic = time.perf_counter()
    axcat.write(results_file)
    toc = time.perf_counter()
    print(f"SAVE in {toc - tic:0.4f} seconds")
    print(axcat)


if __name__ == "__main__":
    help_str = """Path to input json file, File expected to have following format:"""
    parser = argparse.ArgumentParser(description='Integrate w. combat.')
    parser.add_argument('json_file', type=str, help=help_str)
    args = parser.parse_args()
    with open(args.json_file) as f:
        json_data = json.load(f)
    #
    main(json_data)
