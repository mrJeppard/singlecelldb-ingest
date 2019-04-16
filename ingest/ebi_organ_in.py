
import os
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String

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








# Full path to the sqllite db
dbpath = "/home/duncan/work/clusterdb-ingest/cluster.db.ebi.filled.ebicellcounts"
# Path to the data directory filled with ebi anndata objects.
DATADIR = "/home/duncan/work/singlecelldb-ingest/data"

# Connection to the database.
dbstartstr = "sqlite:///%s" % dbpath
engine = create_engine(dbstartstr, echo=False)
Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)
session = Session()

import pandas as pd

def parse_id(filename):
    return filename.split("_")[0]

def make_url(id):
    return "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download?fileType=experiment-design&accessKey=" % id

def get_dataset(name, session):
    return session.query(Dataset.id).filter(Dataset.name == name)


def print_uniques(df, key, c):
    try:
        print(key, df[key].unique().tolist())
    except KeyError:
        c+=1
    return c

for filename in os.listdir(DATADIR):


    df = pd.read_csv(make_url(parse_id(filename)), sep="\t")

    try:
        organs = df['Sample Characteristic[organism part]'].unique().tolist()
        print(organs, filename)
        dataset_id = get_dataset(parse_id(filename), session).update(
            {"organ": ", ".join(organs)}
        )
        session.commit()

    except KeyError:
        print(filename, " did not have organism part with your search")


session.close()