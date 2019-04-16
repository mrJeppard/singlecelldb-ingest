"""
Take in a scanpy object with a cluster solution and write out the ccentroids to csv
"""

"""
Assumes gene names are originally in ensemble space.

"""

import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import sys


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




def centroids2(ad, cs_name="res.1"):
    cluster_solution = ad.obs[cs_name]
    # Calculate each centroid.
    centroids = pd.DataFrame(index=ad.var_names)
    for cluster_name in cluster_solution.unique():
        centroid = pd.Series(ad[ad.obs.index[ad.obs[cs_name] == cluster_name]].X.mean(axis=0).tolist()[0], index=ad.var_names)
        centroids[cluster_name] = centroid

    return centroids

def main():

    ad_path, cs_name, output = parse_args()

    ad = sc.read(ad_path)
    centroids2(ad, cs_name).to_csv(output, header=True)


if __name__ == "__main__":
    sys.exit(main())
