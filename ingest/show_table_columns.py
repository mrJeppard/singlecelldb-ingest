
import sqlite3
import pprint
from sqlalchemy import create_engine, inspect

def show(database_path, table_name="cluster"):
    # Connection to the database.
    dbstartstr = "sqlite:///%s" % database_path
    engine = create_engine(dbstartstr, echo=False)
    inspector = inspect(engine)

    # Get column information
    pprint.pprint(inspector.get_columns('dataset'))

show("/Users/swat/dev/cdbIngest/cluster.db", "dataset")
