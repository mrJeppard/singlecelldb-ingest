
import sqlite3
import pprint
from sqlalchemy import create_engine, inspect, MetaData, select, Table

def show(database_path, table_name="dataset", column="publication_url"):
    # Connection to the database.
    dbstartstr = "sqlite:///%s" % database_path
    engine = create_engine(dbstartstr, echo=False)
    inspector = inspect(engine)
    conn = engine.connect()
    meta = MetaData(engine,reflect=True)
    table = meta.tables[table_name]

    select_st = select([
        table.c.publication_url]).where(
        table.c.publication_url != None)
    res = conn.execute(select_st)
    for _row in res:
        print(_row)


show("/Users/swat/dev/cdbIngest/cluster.db", "dataset", "publication_url")
