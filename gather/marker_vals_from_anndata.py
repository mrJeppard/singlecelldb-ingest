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
from scipy.stats import hypergeom



def hypergeo_z(actual, mu, std):
    return (actual - mu) / std

def hypergeo_mu_std(n_smaller_set, n_larger_set, total_number_of_genes=30000):
    [M, n, N] = [total_number_of_genes, n_smaller_set, n_larger_set]
    hg = hypergeom(M, n, N)
    return hg.mean(), hg.std()

def hypergeo(gene_intersection, n_smaller_set, n_larger_set,  total_number_of_genes=30000):
    return hypergeom.sf(gene_intersection, total_number_of_genes, n_smaller_set, n_larger_set)

def two_averages_fold_change(centroids1, centroids2):
    # subsetted down to have the same features.


    markers1 = find_markers_clusters(centroids1)
    markers2 = find_markers_clusters(centroids2)

    markers1.shape
    markers2.shape

    n_gene_intersection = len(set(markers1.index).intersection(markers2.index))
    gene_union = set(markers1.index).union(markers2.index)
    multiplier = hypergeo_z(n_gene_intersection, *hypergeo_mu_std(len(markers1), len(markers2)))

    centroids1 = centroids[centroids.columns[0:2]]
    centroids2 = centroids[centroids.columns[-2:]]

    second_largest1 = centroids1.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
    second_largest2 = centroids2.apply(lambda x: x.nlargest(2).iloc[1], axis=1)

    largest1 = centroids1.apply(lambda x: x.max(), axis=1)
    largest2 = centroids2.apply(lambda x: x.max(), axis=1)

    largest_cluster1 = centroids1.apply(lambda x: x.idxmax(), axis=1)
    largest_cluster2 = centroids2.apply(lambda x: x.idxmax(), axis=1)

    cluster_pair_ids = pd.Series(
        ["%s:%s" % (f,s) for (f,s) in zip(largest_cluster1, largest_cluster2)],
        index=largest_cluster2.index
    )

    second_largest_combo = second_largest1 + second_largest2
    largest_combo = largest1 + largest2
    cluster_name_df = largest_combo
    second_largest_combo = np.log2(second_largest_combo + 2)
    log_centroids = np.log2(largest_combo + 2)
    two_fc = log_centroids - second_largest_combo

    new_genes = two_fc.index[(two_fc > 1).tolist()]
    new_cluster_pairs = cluster_pair_ids.loc[new_genes]
    cluster_pair_scores = {}
    for cluster_pair_id in new_cluster_pairs.unique():
        cluster1, cluster2 = cluster_pair_id.split(":")
        genes_in_pair = new_cluster_pairs.index[new_cluster_pairs == cluster_pair_id]
        print(len(genes_in_pair))

        try:
            m1 = set(markers1.index[markers1 == cluster1])
            m2 = set(markers1.index[markers2 == cluster2])
            n_gene_intersection = len(m1.intersection(m2))
            multiplier = hypergeo_z(n_gene_intersection, *hypergeo_mu_std(len(m1), len(m2)))
            print(multiplier)
            multiplier = np.max(multiplier, 1)
        except IndexError:
            multiplier = 1

        each_genes = set(markers1.index[markers1 == cluster1]).union(markers2.index[markers2 == cluster2])
        jac = len(each_genes.intersection(genes_in_pair)) / len(each_genes)
        cluster_pair_scores[cluster_pair_id] = jac * multiplier

    return cluster_pair_scores


def log2_fold_change(A,B):
    return math.log(A, 2) - math.log(B, 2)

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-ia',"--anndata", type=str, help="anndata object from scanpy",
                        required=True
                        )

    parser.add_argument('-ic', "--csv", type=str, help="expression counts .csv if different than clustered.ad",
                        required=False, default=None
                        )

    parser.add_argument('-cs', "--cluster_string", type=str, help="anndata.uns key that has a list of all available clusters",
                        required=False
                        )

    parser.add_argument('-o',"--out_dir", type=str, help="name of directory to put gzipped json file",
                        required=True
                        )

    opts = parser.parse_args()

    anndata, counts, cluster_string, out_file = opts.anndata, opts.csv, opts.cluster_string, opts.out_dir

    return anndata, counts, cluster_string, out_file


def write_data(out_file, marker_dicts):
    with gzip.GzipFile(out_file, 'w') as fout:
        fout.write(json.dumps({"markers": marker_dicts}).encode('utf-8'))


def get_counts(ad, counts):
    if isinstance(counts, pd.DataFrame):
        return counts
    return ad.to_df()

def datasetname(filename):
    day = filename.split("-")[0]
    return "hPSC-cardiac-dev_%s" % day


def find_markers_clusters(centroids):
    # Find the second largest of each centroid
    second_largest = centroids.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
    largest = centroids.max(axis=1)
    largest_cluster_name = centroids.apply(lambda x: x.idxmax(), axis=1)
    # Do log transforms with pseudocount
    second_largest = np.log2(second_largest + 2)
    log_centroids = np.log2(largest + 2)
    # Find the log2 fold change from next highest.
    sdiff = log_centroids.sub(second_largest, axis=0)
    # Keep only genes that have at least 1 log fold change (fold change of 2) higher than next.
    marker_genes = sdiff.index[(sdiff >= 1).tolist()]
    cluster_idd = largest_cluster_name[marker_genes]

    return cluster_idd


def find_markers(centroids):
    # Find the second largest of each centroid
    second_largest = centroids.apply(lambda x: x.nlargest(2).iloc[1], axis=1)
    # Do log transforms with pseudocount
    second_largest = np.log2(second_largest + 2)
    log_centroids = np.log2(centroids + 2)
    # Find the log2 fold change from next highest.
    sdiff = log_centroids.sub(second_largest, axis=0)
    # Keep only genes that have at least 1 log fold change (fold change of 2) higher than next.
    marker_genes = sdiff.index[(sdiff.max(axis=1) >= 1).tolist()]
    #print("found %d possible marker genes for cluster solution %s in datasest %s" % (
    #len(marker_genes), cluster_solution_name, dataset_name))

    return marker_genes


def detail_marker_dicts(marker_genes, centroids, n_xpr, n_not_xpr):
    marker_dicts = []
    for gene in marker_genes:
        second_largest = centroids.loc[gene].nlargest(2).tolist()[1] + 2
        minimum = centroids.loc[gene].min() + 2

        for cluster_name in centroids.columns:
            gene_centroid = centroids.loc[gene, cluster_name] + 2
            other_cluster_names = [name for name in centroids.columns if name != cluster_name]
            # TP: true positives
            expressed_in_cluster = n_xpr.loc[gene, cluster_name]
            # FN: false negatives
            not_expressed_in_cluster = n_not_xpr.loc[gene, cluster_name]
            # FP: false positives
            expressed_out_cluster = n_xpr.loc[gene, other_cluster_names].sum()
            # TN: true negatives
            not_expressed_out_cluster = n_not_xpr.loc[gene, other_cluster_names].sum()
            # FP + TN
            out_size = not_expressed_out_cluster + expressed_out_cluster
            # TP + FN
            cluster_size = expressed_in_cluster + not_expressed_in_cluster

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

    return marker_dicts


def main():

    in_file, counts_file, cluster_string, out_dir = parse_args()

    ad = sc.read(in_file)

    try:
        counts = pd.read_csv(counts_file, index_col=0).transpose()
    except ValueError:
        pass

    marker_dicts = []
    dataset_name = datasetname(os.path.basename(in_file))

    out_file = os.path.join(out_dir, "%s_markers.json.gzip" % dataset_name)

    if cluster_string is None:
        cluster_solution_names = ["louvain"]
    else:
        cluster_solution_names = ad.uns[cluster_string]

    for cluster_solution_name in cluster_solution_names:
        X = get_counts(ad, counts)
        X = (X / X.sum()) * 10000
        cluster_solution = ad.obs[cluster_solution_name]
        cluster_solution = cluster_solution.dropna()

        # Calculate each centroid.
        centroids = pd.DataFrame(index=X.columns)
        for cluster_name in cluster_solution.unique():
            cell_names = cluster_solution.index[(cluster_solution == cluster_name).tolist()]
            centroid = X.loc[cell_names].mean(axis=0)
            centroids[cluster_name] = centroid

        marker_genes = find_markers(centroids)
        print("found %d possible marker genes for cluster solution %s in datasest %s" % (len(marker_genes), cluster_solution_name, dataset_name))
        pd.Series(marker_genes).to_csv("trash.txt",index=False)
        # Subset down to only those marker genes.
        X = X[marker_genes]

        # Now go through all the found marker genes and calculate each metric for each cluster.
        print(marker_genes)
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
