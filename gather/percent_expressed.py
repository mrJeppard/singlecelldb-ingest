"""
Outputs a .csv where columns are cluster name and rows are genes with values of the percent of cells that had one
of the marker genes expressed.

"""
import argparse
import pandas as pd
import scanpy as sc
import sys
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-ia',"--anndata", type=str, help="path to anndata object",
                        required=True
                        )

    parser.add_argument('-cs', "--cluster_solution_name", type=str, help="name of the accessor to the cluster solution",
                        required=True
                        )

    parser.add_argument('-o', "--output", type=str, help="name of the centroids .csv file ",
                        required=True
                        )

    opts = parser.parse_args()
    ad_path, cs_name, output = opts.anndata, opts.cluster_solution_name, opts.output

    return ad_path, cs_name, output


def proportion_expressed_cluster(ad, cluster_solution_name):
    """
    outputs a dataframe [genes x cluster] that is the percentage that gene is expressed in a cluster
    :param ad: scanpy.Anndata
    :param cluster_solution_name: string, key accessor for ad.obs cluster_solution (
    :return: pandas.DataFrame
    """

    cluster_solution = ad.obs[cluster_solution_name]

    pcent_df = pd.DataFrame(index=ad.var_names)
    for cluster_name in cluster_solution.unique():
        gt0 = (ad[ad.obs_names[cluster_solution == cluster_name]].X.toarray() > 0).sum(axis=0)
        pcent_df[cluster_name] = gt0 / len(ad.obs_names)

    return pcent_df


def marker_score(proportion_in, proportion_out, pseudocount=.1):
    difference_in_proportions = proportion_in - proportion_out
    scores = np.exp2(difference_in_proportions) / (np.abs(difference_in_proportions) + pseudocount)

    return scores


def proportion_out(df, cluster_name):
    """
    Returns the proportion of cells expressed in all other groupings other than cluster_name.
    :param df:
    :param cluster_name:
    :return:
    """
    return df[[cn for cn in df if cn != cluster_name]].sum(axis=1)


def marker_scores(proportions, pseudocount=.1):
    scores = proportions.apply(lambda x: marker_score(x, proportion_out(proportions, x.name), pseudocount), axis=0)
    return scores


def main():

    ad_path, cs_name, output = parse_args()
    ad = sc.read(ad_path)
    proportions = proportion_expressed_cluster(ad, cs_name)
    marker_scores(proportions).to_csv(output, header=True)


if __name__ == "__main__":
    sys.exit(main())
