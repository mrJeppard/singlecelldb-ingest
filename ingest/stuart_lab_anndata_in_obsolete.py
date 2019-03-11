# add data to the database from an annData scanpy object
import scanpy.api as sc
import pandas as pd
import sklearn



def read_scanpy_object():
    # read in the scanpy object (clustered object and raw/unfiltered gene expression data), for example:
    sc_obj = sc.read("/projects/sysbio/users/cellAtlas/data/cluster/10xGenomics_pbmc8k_clustered.h5ad")
    sc_obj.raw = sc.read("/projects/sysbio/users/cellAtlas/data/scanpyObj/logarithmized_raw_data/10xGenomics_pbmc8k_raw_log.h5ad")


def expression_matrix():
    # write the raw (unfiltered) expression matrix to file

    # if the expression matrix (sc_obj.raw.X) is stored in a sparse format
    expression_df = pd.DataFrame(data=sc_obj.raw.X.toarray(), index = sc_obj.obs.index.tolist(), columns = sc_obj.raw.var.index.tolist())
    expression_df.to_csv("./cluster_db_input_test/pbmc8k_raw_log_expr_matrix.tsv", sep='\t', index=True,header=True)

    # if the expression matrix (sc_obj.raw.X) is stored in a matrix format
    #expression_df = pd.DataFrame(data=sc_obj.raw.X, index = sc_obj.obs.index.tolist(), columns = sc_obj.raw.var.index.tolist())
    #expression_df.to_csv("./cluster_db_input/XXX_raw_log_expr_matrix.tsv", sep='\t', index=True,header=True)


def cell_assignment():
    # write the cell assignments to file
    louvain_cluster_df = pd.DataFrame(data=0, index = sc_obj.obs.louvain.index.tolist(), columns = ["sample","cluster"])
    louvain_cluster_df["sample"] = sc_obj.obs.louvain.index.tolist()
    louvain_clu ster_df["cluster"] = sc_obj.obs.louvain
    louvain_cluster_df.to_csv("./cluster_db_input_test/pbmc8k_louvain100pcs_clustering_table.tsv", sep='\t', index=False,header=True)


def compute_average_gene_expression_per_cluster():
    # compute the average gene expression per cluster
    sorted_k = sorted([int(k) for k in list(set(sc_obj.obs.louvain))])
    avg_expr_df = pd.DataFrame(data=0, index = sc_obj.raw.var.index.tolist(), columns = sorted_k)
    for cluster in sorted_k:
        cluster_specific_cells = sc_obj.obs.loc[(sc_obj.obs.louvain == str(cluster))].index.tolist()
        sample_cluster_subset = sc_obj[cluster_specific_cells,:]
        subset_df = pd.DataFrame(data=sample_cluster_subset.raw.X.toarray(), index = sample_cluster_subset.obs.index.tolist(), columns = sample_cluster_subset.raw.var.index.tolist())
        avg_over_cells = [x/subset_df.shape[0] for x in subset_df.sum(axis=0)]
        avg_expr_df[cluster] = avg_over_cells


def average_gene_expression_per_cluster():
    # write the average gene expression per cluster to file
    avg_expr_df.T.to_csv("./cluster_db_input_test/pbmc8k_louvain100pcs_avg_expr_per_cluster.tsv", sep='\t', index=True,header=True)


def compute_gene_of_set()
    # compute the signature genes (AUC method)
    louvain_clusters = sorted(set(sc_obj.obs.louvain))
    auc_df_louvain = pd.DataFrame(data=0, index = sc_obj.raw.var.index.tolist(), columns = louvain_clusters)
    for c in louvain_clusters:
        print(c)
        cluster_vec = [1 if x==c else 0 for x in sc_obj.obs.louvain]
        full_matrix = sc_obj.raw.X.toarray()
        for gene in sc_obj.raw.var.index.tolist():
            gene_vec = full_matrix[:,list(sc_obj.raw.var.index).index(gene)]
            auc = sklearn.metrics.roc_auc_score(cluster_vec, gene_vec)
            auc_df_louvain.loc[gene,c] = auc

    auc_df_louvain_inverted = auc_df_louvain.copy()
    for c in louvain_clusters:
        auc_df_louvain_inverted[c] = [1-x if x<0.5 else x for x in auc_df_louvain_inverted[c]]

    auc_df_louvain_inverted['max'] = auc_df_louvain_inverted.max(axis=1)
    AUCsig_louvainClusters_genes = list(auc_df_louvain_inverted[auc_df_louvain_inverted['max']>0.75].index)


def gene_of_set():
    # bring signature genes in correct format and write to file
    first = ["gene_set_name","method","clustering_solution_name","dataset_name"] + ["gene_"+str(i) for i in range(len(AUCsig_louvainClusters_genes))]
    second = ["AUCge75","gene included if ROCAUC >= 0.75 for any cluster","louvain100pcs","pbmc8k"] + AUCsig_louvainClusters_genes

    gene_set_df = pd.DataFrame(data=0, index = first, columns = ["second"])
    gene_set_df["second"] = second

    gene_set_df.to_csv("./cluster_db_input_test/pbmc8k_louvain100pcs_AUCge75_sigGeneSet.tsv", sep='\t', index=True,header=False)
