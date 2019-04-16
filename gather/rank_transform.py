"""

"""

"""
Take in a centroids.csv and output a nes_score.csv, input a .gmt file.
"""

"""
Assumes gene names are originally in ensemble space.

"""

import argparse
import gseapy as gs
import numpy as np
import pandas as pd
import sys


def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--centroids", type=str, help="path to the centroids.csv file",
                        required=True
                        )


    parser.add_argument('-o', "--output", type=str, help="name of the nes score .csv file",
                        required=True
                        )

    opts = parser.parse_args()
    centroids_path, output = opts.centroids, opts.output

    return centroids_path, output

def fc_vs_all(centroids):
    fc = pd.DataFrame(index=centroids.index)
    for cname in centroids:
        other_names = [cn for cn in centroids if cn != cname]
        #print(other_names)
        centroid = centroids[cname]
        others_centroid = centroids[other_names]
        #print(others_centroid.head())
        others_centroid = others_centroid.sum(axis=1)
        #print(others_centroid)
        others_centroid = np.log2(others_centroid + 2)
        log_centroid = np.log2(centroid + 2)
        fc[cname] = log_centroid.sub(others_centroid, axis=0)

    #print(fc.head())
    return fc


def fc_vs_next(centroids):
    """
    for each cluster compares the fold change of that gene from the next highest expressing cluster.
    :param centroids: pandas.DataFrame
    :return: pandas.DataFrame (columns clusters rows gene)
    """
    second_largest = centroids.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
    second_largest = np.log2(second_largest + 2)

    fc = pd.DataFrame(index=centroids.index)
    for cname in centroids:
        centroid = centroids[cname]
        log_centroid = np.log2(centroid + 2)
        fc[cname]=log_centroid.sub(second_largest, axis=0)

    return fc

def all_gene_sets(gmt_path):
    gene_sets = []
    with open(gmt_path, "r") as gmt:
        for line in gmt.readlines():
            gene_sets += [line.split("\t")[0]]

    return gene_sets

def main():

    centroids_path, output = parse_args()

    centroids = pd.read_csv(centroids_path, index_col=0)

    ranks_per_cluster = fc_vs_all(centroids)

    ranks_per_cluster[ranks_per_cluster.columns[0]]\
        .to_csv(output, sep="\t", header=True)

if __name__ == "__main__":
    sys.exit(main())
