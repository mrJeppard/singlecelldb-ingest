"""
Clear the database. Good for testing.
"""
import os
from sqlalchemy import create_engine, Table, MetaData, select

# Full path to the sqllite db
FULLDBPATH = "/Users/swat/dev/cdb/cluster.db"

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

# !!!!!!!!!!!! clear the DB for testing.
conn.execute(dataset.delete())
conn.execute(cluster_solution_table.delete())
conn.execute(cluster.delete())
conn.execute(cell_of_cluster.delete())
