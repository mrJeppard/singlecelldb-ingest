"""
Take in a centroids.csv and output a nes_score.csv, input a .gmt file.
"""

"""
Assumes gene names are originally in ensemble space.

"""

import argparse
import gseapy as gp
import numpy as np
import pandas as pd
import sys


def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--centroids", type=str, help="path to the centroids.csv file",
                        required=True
                        )

    parser.add_argument('-g', "--gmt_file", type=str, help="path to the .gmt file",
                        required=True
                        )

    parser.add_argument('-o', "--output", type=str, help="name of the nes score .csv file",
                        required=True
                        )


    opts = parser.parse_args()
    centroids_path, gmt_path, output = opts.centroids, opts.gmt_file, opts.output

    return centroids_path, gmt_path, output

def fc_vs_next(centroids):
    """
    for each cluster compares the fold change of that gene from the next highest expressing cluster.
    :param centroids: pandas.DataFrame
    :return: pandas.DataFrame (columns clusters rows gene)
    """
    genes_considering = centroids.std(axis=1).sort_values().index[-5000:]
    centroids = centroids.loc[genes_considering]

    second_largest = centroids.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
    second_largest = np.log2(second_largest + 2)

    fc = pd.DataFrame(index=centroids.index)
    for cname in centroids:
        centroid = centroids[cname]
        log_centroid = np.log2(centroid + 2)
        fc[cname]=log_centroid.sub(second_largest, axis=0)

    return fc


def fc_vs_all(centroids, psuedocount=.001):
    genes_considering = centroids.std(axis=1).sort_values().index[-5000:]
    centroids = centroids.loc[genes_considering]
    fc = pd.DataFrame(index=centroids.index)

    for cname in centroids:
        other_names = [cn for cn in centroids if cn != cname]
        #print(other_names)
        centroid = centroids[cname]
        others_centroid = centroids[other_names]
        #print(others_centroid.head())
        others_centroid = others_centroid.sum(axis=1)
        #print(others_centroid)
        others_centroid = np.log2(others_centroid + psuedocount)
        log_centroid = np.log2(centroid + psuedocount)
        fc[cname] = log_centroid.sub(others_centroid, axis=0)

    #print(fc.head())
    return fc


def all_gene_sets(gmt_path):
    gene_sets = []
    with open(gmt_path, "r") as gmt:
        for line in gmt.readlines():
            gene_sets += [line.split("\t")[0]]

    return gene_sets

def nes_vectors(ranks, gene_set_file, gene_sets):

    import gseapy as gp
    clusters = ranks.columns.tolist()

    # Just to get pathway names
    nes_df = pd.DataFrame(index=gene_sets)

    for cluster in clusters:
        print("computing for cluster %s" % cluster)
        rank = ranks[[cluster]].reset_index()
        pre_res = gp.prerank(rnk=rank, gene_sets=gene_set_file,
                             processes=10,
                             permutation_num=100,  # reduce number to speed up test
                             outdir="/home/duncan/work/clusterdb-ingest/test/trash")
        nes_df[cluster] = pre_res.res2d["nes"]
        print("number of positive nes scores %d" % (nes_df[cluster] > 0).sum())

    return nes_df

def adj_enrich_score(ranks, gene_set_file, gene_sets, cutoff=.5):
    clusters = ranks.columns.tolist()

    # Just to get pathway names
    nes_df = pd.DataFrame(index=gene_sets)

    for cluster in clusters:
        print("computing for cluster %s" % cluster)
        genes = ranks.index[(ranks[cluster]>cutoff).tolist()]
        #print(genes)
        print(len(genes))

        enr = gp.enrichr(gene_list=genes.tolist(),
                         # or gene_list=glist
                         description='test_name',
                         gene_sets=gene_set_file,
                         outdir='../test/enrichr_kegg',
                         cutoff=1  # test dataset, use lower value of range(0,1)
                         )
        try:
            enr.results.index = enr.results["Term"]
            #print(enr.results.head())
            #print(nes_df.head())
            nes_df[cluster] = enr.results["Adjusted P-value"]
        except KeyError:
            print(enr.results.columns)
            #print(nes_df.head())

        #print("number of positive nes scores %d" % (nes_df[cluster] > 0).sum())

    return nes_df


def main():

    centroids_path, gmt_path, output = parse_args()

    centroids = pd.read_csv(centroids_path, index_col=0)

    ranks_per_cluster = fc_vs_all(centroids)

    print((ranks_per_cluster>0).sum())
    nes_vectors(ranks_per_cluster, gmt_path, gene_sets=all_gene_sets(gmt_path))\
        .to_csv(output, header=True)

if __name__ == "__main__":
    sys.exit(main())
