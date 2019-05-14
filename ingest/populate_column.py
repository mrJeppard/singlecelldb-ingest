
import sqlite3
import pprint
from sqlalchemy import create_engine, inspect, MetaData, select, Table, update

def populate_column(database_path, table_name="dataset"):

    # Connection to the database.
    dbstartstr = "sqlite:///%s" % database_path
    engine = create_engine(dbstartstr, echo=False)
    inspector = inspect(engine)
    conn = engine.connect()
    meta = MetaData(engine,reflect=True)
    table = meta.tables[table_name]

    # update the column in every row.
    query = update(table).values(role='public')
    res = conn.execute(query)
    
    # verify update.
    select_st = select([table.c.role])
    res = conn.execute(select_st)
    for _row in res:
        print(_row)

populate_column('/Users/swat/dev/cdb/cluster.db')
