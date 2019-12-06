"""

"""
import sqlite3

import PeakSQL.tables as tables


class DataBase(object):
    def __init__(self, db='PeakSQL.sqlite'):
        self.db = db

        self.conn = sqlite3.connect(db)
        self.cursor = self.conn.cursor()

        # register all the tables (Species, Chromosome, Condition, Peak)
        for table in [table for table in dir(tables) if not table.startswith('__')]:
            self.cursor.execute(f"CREATE TABLE IF NOT EXISTS {getattr(tables, table)}")

        self.conn.commit()

    # imported methods
    from ._database_sql import add_data, add_assembly, get_id

    @property
    def assemblies(self):
        """
        return all registred assemblies
        """
        return [val[0] for val in self.cursor.execute(f"SELECT Assembly FROM Species").fetchall()]
