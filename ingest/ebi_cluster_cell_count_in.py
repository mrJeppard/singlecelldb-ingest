"""
Works with a directory of anndata objects which are the result of anndataize_ebi.py

To run change the DATADIR and FULLPATHDB global run and run the script from the repos venv.

python ebi_anndata_in.py


"""

import os
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, DECIMAL, ForeignKey, Float

import scanpy as sc

Base = declarative_base()


class Dataset(Base):

    __tablename__ = "dataset"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, nullable=False, unique=True)
    uuid = Column(String)
    species = Column(String)
    organ = Column(String)
    cell_count = Column(Integer)
    disease = Column(String)
    platform = Column(String)
    description = Column(String)
    data_source_url = Column(String)
    publication_url = Column(String)

    def __repr__(self):
     return "<User(name=%s, species=%s, organ=%s, cell count=%d, data source url=%s )>" % \
            (self.name, self.species, self.organ, self.cell_count, self.data_source_url)


class ClusterSolution(Base):

    __tablename__ = "cluster_solution"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, nullable=False)
    description = Column(String)
    method = Column(String)
    method_implementation = Column(String)
    method_url = Column(String)
    method_parameters = Column(String)
    scores = Column(String)
    analyst  = Column(String)
    likes = Column(Integer)
    expression_hash = Column(String)
    dataset_id = Column(Integer, ForeignKey("dataset.id"), nullable=False)


class Cluster(Base):

    __tablename__ = "cluster"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, nullable=False)
    label = Column(String)
    description = Column(String)
    cell_count = Column(Integer)
    cluster_solution_id = Column(Integer, ForeignKey("cluster_solution.id"), nullable=False)


class CellAssignment(Base):

	__tablename__ = "cell_of_cluster"

	id = Column(String, primary_key=True, autoincrement=True)
	name = Column(String, nullable=False)
	cluster = Column(Integer, ForeignKey("cluster.id"), nullable=False)


class Marker(Base):

	__tablename__ = "marker"

	id = Column(Integer, primary_key=True)
	hugo_name = Column(String)
	ensembl_name = Column(String)
	sensitivity = Column(DECIMAL)
	specificity = Column(DECIMAL)
	precision = Column(DECIMAL)
	recall = Column(DECIMAL)
	accuracy = Column(DECIMAL)
	t_pval = Column(DECIMAL)
	z_pval = Column(DECIMAL)
	t_stat = Column(Float)
	z_stat = Column(Float)
	log2_fold_change_vs_min = Column(Float)
	log2_fold_change_vs_next = Column(Float)
	mean_expression = Column(Float)
	dataset_id = Column(Integer, ForeignKey("dataset.id"), nullable=False)
	cluster_solution_id = Column(Integer, ForeignKey("cluster_solution.id"), nullable=False)
	cluster_id = Column(Integer, ForeignKey("cluster.id"), nullable=False)

	def __repr__(self):

		return "<User(name=%s, tstat=%f, mean=%f, log2fc_next=%f)>" % \
			   (self.hugo_name, self.t_stat, self.mean_expression, self.log2_fold_change_vs_next)



# Full path to the sqllite db
dbpath = "/home/duncan/work/clusterdb-ingest/cluster.db.ebi.filled"
# Path to the data directory filled with ebi anndata objects.
DATADIR = "/home/duncan/work/singlecelldb-ingest/data"
# String to access all the clusters vars
CLUSTERSTR = "sc3_cluster_solutions"

# Connection to the database.
dbstartstr = "sqlite:///%s" % dbpath
engine = create_engine(dbstartstr, echo=False)
Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)
session = Session()

def get_dataset_id(name, session):
    return session.query(Dataset.id).filter(Dataset.name == name).one()[0]


def get_cluster_solution_id(dataset_id, name, session):
    return session.query(ClusterSolution.id).filter(
                and_(
                    ClusterSolution.dataset_id == dataset_id,
                    ClusterSolution.name == name
                )
            ).one()[0]


def get_cluster_id(cluster_solution_id, name, session):
    return session.query(Cluster.id).filter(
        and_(
            Cluster.cluster_solution_id == cluster_solution_id,
            Cluster.name == name
        )
    ).one()[0]


def update_with_cell_count(cluster_id, cell_count, session):
    session.query(Cluster).filter(Cluster.id == cluster_id).update({"cell_count": cell_count})
    try:
        session.commit()
    except:
        session.rollback()
        raise



for filename in os.listdir(DATADIR):
    print(filename)
    ad = sc.read(os.path.join(DATADIR,filename))

    name = filename.split("_ebi")[0]
    cell_count = ad.shape[0]
    data_source_url = ad.uns["view_data_url"]
    try:
        cluster_solution_names = ad.uns[CLUSTERSTR]
    except KeyError:
        cluster_solution_names = "louvain"


    cluster_solutions = ad.obs[cluster_solution_names]
    dataset_id = get_dataset_id(name, session)

    for cluster_solution_name in cluster_solution_names:
        cluster_solution_id = get_cluster_solution_id(dataset_id, cluster_solution_name, session)
        cluster_solution = cluster_solutions[cluster_solution_name].dropna()
        cell_counts = cluster_solution.value_counts()
        for cluster_name in cell_counts.index:
            cluster_id = get_cluster_id(cluster_solution_id, cluster_name, session)
            update_with_cell_count(cluster_id, cell_counts[cluster_name].item(), session)


session.close()