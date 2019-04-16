"""
Weighted correlation of NES vectors. The idea is to weight positive values on either vector. Positive values
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def corr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))

def weighted_corr(x, y):
    w = weights(x,y)
    return corr(x,y,w)


def all_by_all(df, weights=None):
    if weights is None:
        weights = df.std(axis=1).values.transpose()

    sim_func = lambda x, y: corr(x, y, weights)
    wd = pdist(df.transpose().values, metric=sim_func)
    wd = squareform(wd)
    wd = pd.DataFrame(wd, columns=df.columns, index=df.columns)
    return wd



def std_weights(invivo_nes, invitro_nes):
    std1 = invivo_nes.std(axis=1)
    std2 = invitro_nes.std(axis=1)
    return np.sqrt(std1.multiply(std2))




invivo_nes = pd.read_csv("./tmp/invivo.celltype.nes.csv", index_col=0)
#invitro_nes = pd.read_csv("./friedmanC.celltype.nesscores.csv", index_col=0)
invitro_nes = pd.read_csv("./tmp/invitro.res_0_4.nes.csv", index_col=0)
#invitro_nes = pd.read_csv("./tmp/test2.tab", index_col=0)

invivo_nes.head()
invivo_nes.shape
invivo_nes.isna().sum()

invitro_nes.shape
invitro_nes.isna().sum()

invivo_nes = invivo_nes.dropna(axis=0, how="any")
invivo_nes.shape

invitro_nes = invitro_nes.dropna(axis=0, how="any")
invitro_nes.shape

pathway_intersect = list(set(invivo_nes.index).intersection(invitro_nes.index))
len(pathway_intersect)
invivo_nes = invivo_nes.loc[pathway_intersect]
invitro_nes = invitro_nes.loc[pathway_intersect]

nes = pd.concat([invivo_nes, invitro_nes], axis=1)
nes.head()

corrs = all_by_all(nes, std_weights(invivo_nes, invitro_nes))
corrs.to_csv("./corrs.invivo-celltype.invitro-res-04.csv")
corrs.columns
#corrs["aCM"][invitro_nes.columns].sort_values()
corrs["pluripotent"][invivo_nes.columns].sort_values(ascending=False)
corrs["definitive_cardiomyocyte"][invivo_nes.columns].sort_values(ascending=False)
corrs["non_contractile"][invivo_nes.columns].sort_values(ascending=False)
#corrs["aCM"].sort_values()

################################################
dataset_name = "fetal combined heart of cells"
cluster_solution_name = "heart cell types"
size = "similarity"
color = "MYL7"

color_centroids = pd.read_csv("invitroCombined.res_0_4.centroids.csv", index_col=0).loc[color]
#color_centroids = pd.read_csv("fetalCombined.celltype.centroids.csv", index_col=0).loc[color]
cluster_cell_counts = pd.read_csv("invitroCombined.cluster.cellcounts.csv", index_col=0)

other_dataset = "in vitro combined heart of cells"
other_species = "human"
other_study = "in vitro"
other_organ = "heart"
other_cluster_solution_name = "louvain resolution 0.4"
big_dict = {
    "dataset_name": dataset_name,
    "cluster_solution_name": cluster_solution_name,
    "size_by": size,
    "color_by": color,
    "cluster_similarities": []
}

for celltype in invivo_nes:
    cs_dict = {
        "dataset": {
          "name": other_dataset,
            "species": other_species,
            "organ": other_organ,
            "study": other_study
        },
        "compared_to_cluster": celltype,
        "cluster_solution_name": other_cluster_solution_name,
        "clusters": []

    }
    for cluster in invitro_nes:
        cluster_dict = {
            "name": cluster,
            "size": corrs.loc[celltype, cluster].item(),
            "color": color_centroids[cluster].item(),
            "cell_count": cluster_cell_counts.loc[int(cluster)].item()
        }
        cs_dict["clusters"].append(cluster_dict)

    big_dict["cluster_similarities"].append(cs_dict)
###################

import json
with open('example-similarities.json', 'w') as outfile:
    json.dump(big_dict, outfile, indent=4)

    print(celltype)

#################
import scanpy as sc
filenames_and_days = [
    ("./clusterd-friedman/day0-rep12.h5ad","Day0"),
    ("./clusterd-friedman/day15-rep12.h5ad", "Day15"),
    ("./clusterd-friedman/day2-rep12.h5ad", "Day2"),
    ("./clusterd-friedman/day30-rep12.h5ad", "Day30"),
    ("./clusterd-friedman/day5-rep12.h5ad", "Day5"),
]

dfs = []
for adfile, day in filenames_and_days:
    ad = sc.read(adfile)
    df = ad.obs[["cell_type"]]
    df = df.dropna(axis=1)
    df.index = ["%s-%s" % (name, day) for name in df.index]
    dfs+=[df]

ctypes  = pd.concat(dfs, axis=0)
ctypes.head()
ctypes.to_csv("friedman-celltype-assignments.csv", header=True)
ctypes.head()
len(ctypes.index) == len(ctypes.index.unique())
"""
{
    dataset_name: "fetal combined heart of cells",
	cluster_solution_name: “heart cell types”
	size_by: 'similarity',
	color_by: ‘MYL7',
            cluster_similarities: [
{
    	dataset: {
name: 'in vitro heart combined',
species: 'Homo sapiens',
organ: 'heart',
study: 'in vitro',,
},
compared_to_cluster: “cell type 1”,
    	cluster_solution_name: 'louvain resolution .25',
    	clusters: [
        	{
            	name: 'A',
            	size: 15,
            	color: -0.75,
		cell_count: 34,
        	},
        	{
            	name: 'B',
            	size: 30,
            	color: 0.5,
	cell_count: 88
             	},
	...
    		],
	},
{
	(repeated...)
},
...
]
}



#############################
# Dont go below here....
##############################
invitro_nes2.head()
n1 = invivo_nes1[invivo_nes1.columns[0]]
n2 = invitro_nes2[invitro_nes2.columns[2]]
c = invivo_nes1.corr()
c[c.columns[0]].sort_values()
c['2']
c.head()
len(set(n1.index).intersection(n2.index))


n1 = n1[list(set(n1.index).intersection(n2.index))]
n2 = n2[list(set(n1.index).intersection(n2.index))]
##############################

import numpy as np

# Without weight.
np.corrcoef(n1,n2)


x, y  = n1, n2

(x>0).sum()
(y>0).sum()

positive_ind = np.logical_or(x>0, y>0)

positive_ind.sum()
weights = np.zeros((len(positive_ind),1))
weights[positive_ind] = 10 / positive_ind.sum()
weights[~positive_ind] = 1 / (len(positive_ind) - positive_ind.sum())
# alternative paying attention to stdard devaition across clusters
weights = nes.std(axis=1).values.transpose()
# Gives you way smaller
#
(100 / 5)
weights.sum()

1 - corr(n1.values, n2.values, weights.transpose())


sim_func = lambda x, y: 1 - corr(x,y, weights)
wd = pdist(nes.transpose().values, metric=sim_func)
wd = squareform(wd)
wd.shape
wd
wd = pd.DataFrame(wd, columns=invivo_nes.columns, index=invivo_nes.columns)
wd.head()
wd["aCM"].sort_values()



nes.head()
aa = all_by_all(nes)

aa.head()

aa["aCM"].sort_values()



ab["aCM"][invitro_nes.columns].sort_values()

closest = ab["aCM"][invitro_nes.columns].sort_values().index[0:2]
ad.obs["is"] = np.logical_or(ad.obs["louvain"])
aCMMarkers = ["res.25", "NPPA", "PAM", "MYL7"]
import scanpy as sc
ad = sc.read("in_vitro_clustered.h5ad")
sc.tl.louvain(ad, resolution=.25, key_added="res.25")
sc.pl.umap(ad, color= aCMMarkers)

dataset_name = "in vivo heart of cells"
cluster_solution_name = "heart cell types"
size = "similarity"
color = "NPPA"
other_dataset_name = "in vitro heart of cells"
other_cluster_solution_name = "louvain resolution .25"
import pandas as pd
centroids_invivo = pd.read_csv("fetalCombined.celltype.centroids.csv", index_col=0)
centroids_invitro = pd.read_csv("invitroCombined.res_0_4.centroids.csv", index_col=0)

possible_markers = ["MYL7", "PAM"]
import scanpy as sc
import numpy as np
#ad = sc.read("tm")
len(ad.obs["celltype"].unique())
np.sum(ad[:, possible_markers[0]].X > 0)

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
        fc[cname] = log_centroid.sub(others_centroid, axis=possible_markers)

    #print(fc.head())
    return fc


centroids_invitro.loc[possible_markers]
centroids_invivo.loc[color]
centroids_invivo.index
invivo_fc = fc_vs_all(centroids_invivo)
invitro_fc = fc_vs_all(centroids_invitro)

invitro_fc.loc[possible_markers].max(axis=1)
invivo_fc.loc[color]
"""
{
	cluster_solution_name: “heart cell types”
	size_by: 'similarity',
	color_by: ‘MYL7',
            cluster_similarities: [
{
    	dataset: {
name: 'in vitro heart combined',
species: 'Homo sapiens',
organ: 'heart',
study: 'in vitro',,
},
compared_to_cluster: “cell type 1”,
    	cluster_solution_name: 'louvain resolution .25',
    	clusters: [
        	{
            	name: 'A',
            	size: 15,
            	color: -0.75,
		cell_count: 34,
        	},
        	{
            	name: 'B',
            	size: 30,
            	color: 0.5,
	cell_count: 88
             	},
	...
    		],
	},
{
	(repeated...)
},
...
]
}
"""
"""
{
	gene: 'ALK', / cluster: 'userCluster1',
	size_by: 'sensitivity',
	color_by: 'z_stat',
            cluster_solutions: [
{
    	dataset: {
name: 'dataset name 2',
species: 'Homo sapiens',
organ: 'heart',
study: 'in vivo',,
},
    	cluster_name: 'solution 2',
    	clusters: [
        	{
            	name: 'A',
            	size: 15,
            	color: -0.75,
		cell_count: 34,
        	},
        	{
            	name: 'B',
            	size: 30,
            	color: 0.5,
	cell_count: 88
             	},
	...
    		],
	},
{
	(another dataset/cluster-solution)
},
...
]
}
"""
