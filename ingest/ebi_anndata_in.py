"""
Works with a directory of anndata objects which are the result of anndataize_ebi.py

To run change the DATADIR global and
"""
import os
import scanpy as sc
from sqlalchemy import create_engine, Table, MetaData

# Full path to the sqllite db
FULLDBPATH = "/home/duncan/work/singlecelldb-ingest"
# Path to the data directory filled with ebi anndata objects.
DATADIR = "./data"


# Connection to the database.
dbstartstr = "sqlite:///%s" % FULLDBPATH
engine = create_engine(dbstartstr, echo=True)
metadata = MetaData()
conn = engine.connect()

# Accessor for each of the tables.
dataset = Table('dataset', metadata, autoload=True, autoload_with=engine)
cluster_solution_table = Table('cluster_solution', metadata, autoload=True, autoload_with=engine)
cluster = Table('cluster', metadata, autoload=True, autoload_with=engine)
cell_of_cluster = Table('cell_of_cluster', metadata, autoload=True, autoload_with=engine)


def cluster_description(k):
    return "sc3 clusters k=%d from ebi's single cell atlas" % k


def sc3_method_url():
    return "http://bioconductor.org/packages/release/bioc/html/SC3.html"


for filename in os.listdir(DATADIR):
    print(filename)
    ad = sc.read(os.path.join(DATADIR,filename))

    name = filename.split("_ebi")[0]
    cell_count = ad.shape[0]
    data_source_url = ad.uns["view_data_url"]
    cluster_solution_names = ad.uns["sc3_cluster_solutions"]

    try:
        preferred_cluster_solution = ad.uns["sc3_preferred_cluster"]
    except KeyError:
        preferred_cluster_solution = None

    species = ad.uns["species"]
    description = ad.uns["short_description"]

    cluster_solutions = ad.obs[cluster_solution_names]
    dataset_ins = dataset.insert().values(
            name=name, 
            species=species, 
            cell_count=cell_count,
            description=description,
            data_source_url=data_source_url
            )

    result = conn.execute(dataset_ins)
    dataset_key = result.inserted_primary_key

    for cluster_solution_name in cluster_solution_names:
        cluster_solution = cluster_solutions[cluster_solution_name].dropna()

        cluster_values = cluster_solution.unique().tolist()
        print("**************************************")
        print(cluster_values)
        k = len(cluster_values)
        cluster_sol_ins = cluster_solution_table.insert().values(
                name=cluster_solution_name,
                description=cluster_description(k),
                method="sc3",
                method_url=sc3_method_url(),
                dataset_id=dataset_key[0]
                )

        result = conn.execute(cluster_sol_ins)
        cluster_sol_key = result.inserted_primary_key
        
        for cluster_value in cluster_values:
            cluster_ins = cluster.insert().values(
                    name=str(cluster_value),
                    cluster_solution_id=cluster_sol_key[0]
                    )
            result = conn.execute(cluster_ins)
            cluster_key = result.inserted_primary_key
            cell_ids = cluster_solution[cluster_solution == cluster_value].index
            cells = [dict(name=n, cluster_id=cluster_key[0]) for n in cell_ids]
            print(cells)
            conn.execute(cell_of_cluster.insert(), cells)
