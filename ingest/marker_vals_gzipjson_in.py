"""
Assumes gene names are originally in ensemble space.

"""

import argparse
import gzip
import json
import mygene
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, DECIMAL, ForeignKey, Float
import sys

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-db',"--database", type=str, help="path to sqllite database",
                        required=True
                        )

    parser.add_argument('-gz', "--gzipjson", type=str, help="path to output file of marker_vals_from_anndata",
                        required=True
                        )

    opts = parser.parse_args()
    database, gzipjson = opts.database, opts.gzipjson

    return database, gzipjson


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



def query_datasetid(dataset_name, session):
	return session.query(Dataset.id).filter(Dataset.name == dataset_name).one()[0]


def query_cluster_solution_id(dataset_id, cluster_solution_name, session):
	return session.query(ClusterSolution.id).filter(
		and_(
			ClusterSolution.dataset_id == dataset_id,
			ClusterSolution.name == cluster_solution_name
		)
	).one()[0]


def query_clusterid(cluster_solution_id, cluster_name, session):
	return session.query(Cluster.id).filter(
				and_(
					Cluster.cluster_solution_id == cluster_solution_id,
					Cluster.name == cluster_name
			)
	).one()[0]


def main():

	dbpath, in_file = parse_args()
	# Full path to the sqllite db
	# dbpath = "/home/duncan/work/ebiscanpy-clusterdb/cluster.db.ebi.filled"

	# Connection to the database.
	dbstartstr = "sqlite:///%s" % dbpath
	engine = create_engine(dbstartstr, echo=False)
	Base.metadata.create_all(engine)
	Session = sessionmaker(bind=engine)
	session = Session()

	with gzip.GzipFile(in_file, 'r') as fin:
		data = json.loads(fin.read().decode('utf-8'))

	print("adding %d markers to the database at %s" % (len(data["markers"]), dbpath))
	# Grab mappings for all of your genes and make an ensemble -> symbol dictionary with everything you got hits for.
	mg = mygene.MyGeneInfo()
	# Pull out ensemble gene names
	gene_names = list(set([d["gene_name"] for d in data["markers"]]))
	ginfo = mg.querymany(gene_names, scopes='ensembl.gene', as_dataframe=True)
	ginfo = ginfo.iloc[(ginfo["notfound"] != True).tolist()]
	gene_lookup = dict(zip(ginfo.index, ginfo["symbol"]))

	Markers = []
	for marker_dict in data["markers"]:
		dataset_id = query_datasetid(marker_dict["dataset_name"], session)
		cluster_solution_id = query_cluster_solution_id(dataset_id, marker_dict["cluster_solution_name"], session)
		cluster_id = query_clusterid(cluster_solution_id, marker_dict["cluster_name"], session)

		try:
			hugo_name = gene_lookup[marker_dict["gene_name"]]
		except KeyError:
			hugo_name = None

		Markers += [
			Marker(
				hugo_name = hugo_name,
				ensembl_name = marker_dict["gene_name"],
				sensitivity = marker_dict["sensitivity"],
				specificity = marker_dict["specificity"],
				precision = marker_dict["precision"],
				recall = marker_dict["recall"],
				accuracy = marker_dict["accuracy"],
				t_pval = marker_dict["t_pval"],
				z_pval = marker_dict["z_pval"],
				z_stat = marker_dict["z_stat"],
				t_stat = marker_dict["t_stat"],
				log2_fold_change_vs_next = marker_dict["log2_change_vs_next"],
				log2_fold_change_vs_min = marker_dict["log2_change_vs_min"],
				mean_expression = marker_dict["mean_expression"],
				dataset_id = dataset_id,
				cluster_solution_id = cluster_solution_id,
				cluster_id = cluster_id
			)
		]


	print("adding to database")
	session.close()

	Session = sessionmaker(bind=engine)
	session = Session()

	session.add_all(Markers)
	session.commit()

	session.close()

	print("an example gene is %s" % marker_dict["gene_name"])


if __name__ == "__main__":
	sys.exit(main())
