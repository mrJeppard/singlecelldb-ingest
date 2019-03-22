"""
Input a anndata object and cluster string, output a gzipped json file that is:

{
    marker_genes: [
        {
            "gene_name" : gene,
            "dataset_name": dataset_name,
            "cluster_solution_name": cluster_solution_name,
            "cluster_name": cluster_name,
            "sensitivity": sensivity, #(+1 read considered a guess) / n cells IN cluster
            "specificity": specificity, #(+1 read considered a guess) (n reads < 1 out of cluster / n cells OUT of cluster)
            "precision": precision, #n +1 reads in cluster / n +1 reads out of cluster
            "accuracy": accuracy,
            "recall": recall, #n +1 reads in cluster / n 0 reads in cluster.
            "t_pval": tpval,
            "z_pval": zpval,
            "z_stat": zstat,
            "t_stat": tstat, #(t statistic of gene expression in / out of cluster)
            "log2_change_vs_min": fold_change_min, #( min cluster.... mincluster will always be 0.
            "log2_change_vs_next": fold_change_next, #(next is the level of the second highest cluster, there for the 2nd cluster will always be 1.
            "mean_expression": centroids.loc[gene, cluster_name]
        }, ...
    ]
}

The data set name is currently inferred from the filename, with filename.split("_ebi")[0]

For each of the clustering solutions, any cell with a value of "NA" or "nan" is discarded, and not considered for
the in/out cluster calculations.

The log fold changes use a pseudocount of 2. This controls the log fold change when you have small values. For instance
if the smaller average expression was .01, and the largest was 2, without a pseudocount the fold change would be
large. Adding a pseudocount of 2 means you would be comparing the fold change between 2.01 and 4 instead, which is less
that a single log2 fold change... that is what we want because we don't want to favor small changes in expression.

For the Ztest, a +1 read in a cell is considered a positive indication, the proportion of +1 reads in the cluster
is compared with the proportion of +1 reads out of the cluster with the alternative hypothesis being that the proporiton
in the cluster is larger than outside of the cluster.

The t-statistics are done with the in cluster vs out of cluster.
"""
import argparse
import sys
from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import ttest_ind
import scanpy as sc
import pandas as pd
import math
import numpy as np
import json
import gzip
import os

def log2_fold_change(A,B):
    return math.log(A, 2) - math.log(B, 2)

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-ia',"--anndata", type=str, help="anndata object from scanpy",
                        required=True
                        )

    parser.add_argument('-cs', "--cluster_string", type=str, help="anndata.uns key that has a list of all available clusters",
                        required=True
                        )

    parser.add_argument('-o',"--out_dir", type=str, help="name of directory to put gzipped json file",
                        required=True
                        )

    opts = parser.parse_args()
    anndata, cluster_string, out_file = opts.anndata, opts.cluster_string, opts.out_dir

    return anndata, cluster_string, out_file


def write_data(out_file, marker_dicts):
    with gzip.GzipFile(out_file, 'w') as fout:
        fout.write(json.dumps({"markers": marker_dicts}).encode('utf-8'))


def main():

    in_file, cluster_string, out_dir = parse_args()

    ad = sc.read(in_file)

    marker_dicts = []
    dataset_name = os.path.basename(in_file).split("_ebi")[0]

    out_file = os.path.join(out_dir, "%s_markers.json.gzip" % dataset_name)
    for cluster_solution_name in ad.uns[cluster_string]:
        X = ad.to_df()
        cluster_solution = ad.obs[cluster_solution_name]
        cluster_solution = cluster_solution.dropna()

        # Calculate each centroid.
        centroids = pd.DataFrame(index=ad.var_names)
        for cluster_name in cluster_solution.unique():
            cell_names = cluster_solution.index[(cluster_solution == cluster_name).tolist()]
            centroid = X.loc[cell_names].mean(axis=0)
            centroids[cluster_name] = centroid

        # Find the second largest of each centroid
        second_largest = centroids.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
        # Do log transforms with pseudocount
        second_largest = np.log2(second_largest + 2)
        log_centroids = np.log2(centroids + 2)
        # Find the log2 fold change from next highest.
        sdiff = log_centroids.sub(second_largest, axis=0)
        # Keep only genes that have at least 1 log fold change (fold change of 2) higher than next.
        marker_genes = sdiff.index[(sdiff.max(axis=1) >= 1).tolist()]
        print("found %d possible marker genes for cluster solution %s in datasest %s" % (len(marker_genes), cluster_solution_name, dataset_name))
        
        # Subset down to only those marker genes.
        X = X[marker_genes]

        # Now go through all the found marker genes and calculate each metric for each cluster.
        for gene in marker_genes:
            second_largest = centroids.loc[gene].nlargest(2).tolist()[1] + 2
            minimum = centroids.loc[gene].min() + 2

            for cluster_name in cluster_solution.unique():
                gene_centroid = centroids.loc[gene, cluster_name] + 2
                cell_names = cluster_solution.index[(cluster_solution == cluster_name).tolist()]
                other_cell_names = cluster_solution.index[(cluster_solution != cluster_name).tolist()]

                # TP: true positives
                expressed_in_cluster = (X.loc[cell_names, gene] > 0).sum()
                # FN: false negatives
                not_expressed_in_cluster = (X.loc[cell_names, gene] == 0).sum()
                # FP: false positives
                expressed_out_cluster = (X.loc[other_cell_names, gene] > 0).sum()
                # TN: true negatives
                not_expressed_out_cluster = (X.loc[other_cell_names, gene] > 0).sum()
                # FP + TN
                out_size = len(other_cell_names)
                # TP + FN
                cluster_size = len(cell_names)

                # TP / (FP + TP)
                precision = expressed_in_cluster / (expressed_in_cluster + expressed_out_cluster)
                # TP / (FN + TP)
                recall = expressed_in_cluster / (not_expressed_in_cluster + expressed_in_cluster)

                # TP / ( TP + FN)
                sensivity = expressed_in_cluster / cluster_size
                # TN / (FP + TN)
                specificity = not_expressed_out_cluster / out_size

                # Accuracy = (TP+TN)/(TP+TN+FP+FN)
                accuracy = (expressed_in_cluster + not_expressed_out_cluster) / (out_size + cluster_size)

                fold_change_next = log2_fold_change(gene_centroid, second_largest)
                fold_change_min = log2_fold_change(gene_centroid, minimum)

                zstat, zpval = proportions_ztest(
                    count=[expressed_in_cluster, expressed_out_cluster],
                    nobs=[cluster_size, out_size],
                    alternative='larger'
                )

                tstat, tpval = ttest_ind(X.loc[cell_names, gene], X.loc[other_cell_names, gene])

                marker_dict = {
                    "gene_name": gene,
                    "dataset_name": dataset_name,
                    "cluster_solution_name": cluster_solution_name,
                    "cluster_name": str(cluster_name),
                    "sensitivity": sensivity.item(),  # (+1 read considered a guess) / n cells IN cluster
                    "specificity": specificity.item(),
                    "precision": precision.item(),  # n +1 reads in cluster / n +1 reads out of cluster
                    "accuracy": accuracy.item(),
                    "recall": recall.item(),  # n +1 reads in cluster / n 0 reads in cluster.
                    "t_pval": tpval.item(),
                    "z_pval": zpval.item(),
                    "z_stat": zstat.item(),
                    "t_stat": tstat.item(),  # (t statistic of gene expression in / out of cluster)
                    "log2_change_vs_min": fold_change_min,  # ( min cluster.... mincluster will always be 0.
                    "log2_change_vs_next": fold_change_next,
                    "mean_expression": centroids.loc[gene, cluster_name].item()
                }

                marker_dicts += [marker_dict]

    write_data(out_file, marker_dicts)


if __name__ == "__main__":
    sys.exit(main())