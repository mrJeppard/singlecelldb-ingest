"""
Works with a directory of anndata objects which are the result of the stuart lab
runs from October 2018

Update the data values common across the cluster solution, such as
cluster_description, method...

To run change the DATADIR and FULLPATHDB global run and run the script from the repos venv.

python stuart_lab_anndata_in.py
"""
import os, csv
import scanpy as sc
import pandas as pd
from sqlalchemy import create_engine, Table, MetaData, select

"""
# Full path to the sqllite db on bop
FULLDBPATH = "/soe/swat/cellAtlas/data/cluster.swat.db"
# Path to the data directory filled with anndata objects.
DATADIR = "/soe/swat/cellAtlas/data/cluster"
# Path to the dataset tsv file.
DATASETPATH = "/soe/swat/clusterdb-ingest/dataset.tsv"
"""

# Full path to the sqllite db for swat
FULLDBPATH = "/Users/swat/dev/cdb/cluster.db"
# Path to the data directory filled with anndata objects.
DATADIR = "/Users/swat/dev/cdbIngest/dataIn"
# Path to the dataset tsv file.
DATASETPATH = "/Users/swat/dev/cdbIngest/dataIn/dataset.tsv"

# Connection to the database.
dbstartstr = "sqlite:///%s" % FULLDBPATH
engine = create_engine(dbstartstr, echo=False)
metadata = MetaData()
conn = engine.connect()

# Accessor for each of the tables.
dataset = Table('dataset', metadata, autoload=True, autoload_with=engine)
cluster_solution_table = Table('cluster_solution', metadata, autoload=True, autoload_with=engine)
cluster = Table('cluster', metadata, autoload=True, autoload_with=engine)
cell_of_cluster = Table('cell_of_cluster', metadata, autoload=True, autoload_with=engine)

# Remove obsolete datasets from the db.
with open(DATASETPATH, 'rU') as fin:
    fin = csv.DictReader(fin, delimiter='\t')
    for row in fin:
        dataset_delete = dataset.delete().where(dataset.c.name == row['name'])
        result = conn.execute(dataset_delete)
        print('deleted:', row['name'])
for name in [
    'Tabula Muris facs',
    'Tabula Muris droplet',
    'Quake Brain',
    'UCSC Human Cortex',
    'Immune Bone',
    'Immune Cord']:
        dataset_delete = dataset.delete().where(dataset.c.name == name)
        result = conn.execute(dataset_delete)
        print('deleted:', name)

# Remove obsolete cluster solutions from the db
