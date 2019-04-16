"""
Outputs a .csv where columns are cluster name and rows are genes with values of the average percent of markers expressed
in the cells of a cluster.

"""
import argparse
import pandas as pd
import scanpy.api as sc
import sys
from decorator import decorator

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-ia',"--anndata", type=str, help="path to anndata object",
                        required=True
                        )

    parser.add_argument('-cs', "--cluster_solution_name", type=str, help="name of the accessor to the cluster solution",
                        required=True
                        )

    parser.add_argument('-g', "--gmt_file", type=str, help="path to the .gmt file",
                        required=True
                        )


    parser.add_argument('-o', "--output", type=str, help="name of the centroids .csv file ",
                        required=True
                        )

    opts = parser.parse_args()
    ad_path, cs_name, output, gmt_file = opts.anndata, opts.cluster_solution_name, opts.output, opts.gmt_file

    return ad_path, cs_name, output, gmt_file


def mean_n_expressed_in_set(ad, in_cluster, gene_set):
    """Mean number of markers expressed in a 'in_cluster' set of cells"""
    return (ad[gene_set, in_cluster].X.toarray() > 0).sum(axis=0).mean() / len(gene_set)


def percentage_markers_expressed_in_cluster(ad, cluster_solution_name, gene_dict):
    """
    outputs a dataframe [genes x cluster] that is the percentage that gene is expressed in a cluster
    :param ad: scanpy.Anndata
    :param cluster_solution_name: string, key accessor for ad.obs cluster_solution (
    :return: pandas.DataFrame
    """

    cluster_solution = ad.obs[cluster_solution_name]

    pcent_df = pd.DataFrame(index=gene_dict.keys())
    for cluster_name in cluster_solution.unique():
        in_cluster = ad.obs_names[cluster_solution == cluster_name]
        av = pd.Series(index=gene_dict.keys())

        for gene_set_name, genes in gene_dict.items():
            av[gene_set_name] = mean_n_expressed_in_set(ad, in_cluster, genes)

        pcent_df[cluster_name] = av

    return pcent_df



@decorator
def filename_to_fileobj(func, *args, **kwargs):
    """To use this your file name/ object has to be the first argument of your function."""
    if type(args[0]) == str:
        args = (open(args[0], "r"),) + args[1:]

    return(func(*args, **kwargs))


@filename_to_fileobj
def read_gmt(fileobj):
    """Read a .gmt file into a gene set name -> [genes] dictionary.

    Args:
       fileobj (string or file like object):

    Returns:
       (dict): branch-id -> markers dictionary
    """

    # Currently assumes gmt file.
    gene_sets = {}

    for line in fileobj:
        row = line.strip().split("\t")
        gene_sets[row[0]] = row[2:]

    return gene_sets


@filename_to_fileobj
def read_tsv_gene_set(fileobj):
    """Read a .tsv file containing gene sets name -> [genes] dictionary.

    Args:
       fileobj (string or file like object):

    Returns:
       (dict): branch-id -> markers dictionary
    """

    # Currently assumes gmt file.
    gene_sets = {}

    for line in fileobj:
        row = line.strip().split("\t")
        gene_sets[row[0]] = row[1:]

    return gene_sets

def main():

    ad_path, cs_name, output, gmt_file = parse_args()

    gene_sets = read_gmt(gmt_file)

    ad = sc.read(ad_path)
    percentage_markers_expressed_in_cluster(ad, cs_name, gene_sets).to_csv(output, header=True)


if __name__ == "__main__":
    sys.exit(main())
