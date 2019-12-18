"""

"""
import sqlite3

import PeakSQL.tables as tables
import os

import pybedtools


class DataBase(object):
    def __init__(self, db='PeakSQL.sqlite'):
        self.db = db

        self.conn = sqlite3.connect(db)
        self.cursor = self.conn.cursor()

        self.open_files = {}

        # register all the tables (Assembly, Chromosome, Condition, Peak)
        for table in [table for table in dir(tables) if not table.startswith('__')]:
            self.cursor.execute(f"CREATE TABLE IF NOT EXISTS {getattr(tables, table)}")

        self.conn.commit()

    @property
    def assemblies(self):
        """
        return all registred assemblies
        """
        return [val[0] for val in self.cursor.execute(f"SELECT Assembly FROM Assembly").fetchall()]


    def add_assembly(self, assembly_file, assembly=None, species=None):
        """

        """
        # set defaults if none provided
        assembly = assembly if assembly else os.path.basename(assembly_file).split('.')[0]
        species = species if species else assembly

        # TODO: what should default behaviour be?
        if assembly in self.assemblies:
            return

        # Make sure the assembly hasn't been added yet
        assert assembly not in self.assemblies, f"Assembly '{assembly}' has already been added " \
                                                f"to the database!"

        self.cursor.execute(f"INSERT INTO Assembly (Assembly, Species, Abspath) "
                            f"VALUES ('{assembly}', '{species}', '{assembly_file}')")

        assembly_id = self.get_id('Assembly', 'Assembly', assembly)

        # now load the chromosome sizes
        fasta = pyfaidx.Fasta(assembly_file)
        self.open_files[assembly_file] = fasta
        for sequence_name in fasta.keys():
            sequence = fasta[sequence_name]
            self.cursor.execute(f"INSERT INTO Chromosome (AssemblyId, Size, Chromosome, Sequence) "
                                f"    VALUES('{assembly_id}', '{len(sequence)}', "
                                f"           '{sequence_name}', '{sequence}')")

        self.conn.commit()


    def add_data(self, data_file, assembly, condition=None):
        """

        """
        # check for supported filetype
        *_, extension = os.path.splitext(data_file)
        assert extension in ['.narrowPeak'], f"The file extension you choose is not supported, TODO"

        # check if species it belongs to has already been added to the database
        species_id = self.get_id('Assembly', 'AssemblyId', assembly)
        print(assembly, species_id)
        assert species_id, f"Assembly '{assembly}' has not been added to the database yet."

        # check if the condition does not already exist for this species
        self.cursor.execute(f"SELECT Condition FROM Condition "
                            f"WHERE  SpeciesId='{species_id}' "
                            f"AND    Condition='{condition}'").fetchone()
        # TODO: what should default behaviour be?
        print(self.cursor.fetchone())
        print(self.cursor.fetchone())
        if self.cursor.fetchone() is not None:
            return
        assert not self.cursor.fetchone(), f"Condition '{condition}' has already been added to " \
                                           f"the database for assembly '{assembly}'."

        # add the condition
        self.cursor.execute(f"INSERT INTO Condition VALUES(NULL, '{condition}', '{species_id}')")

        condition_id = self.get_id('Condition', 'Condition', condition)
        bed = pybedtools.BedTool(data_file)
        for region in bed:
            chromosome_id = self.get_id('Chromosome', 'Chromosome', region.chrom)
            if len(region.fields) == 3:
                raise NotImplementedError
            else:
                assert len(region.fields) == 10, "TODO error message"
                region.fields[-1] = int(region.fields[1]) + int(region.fields[-1])
                self.cursor.execute(f"""INSERT INTO Bed """
                                    f"""  VALUES(NULL, '{species_id}', '{condition_id}', """
                                    f"""  '{chromosome_id}', '{"', '".join(region.fields[1:])}')""")

        self.conn.commit()


    def get_id(self, table, column, value):
        """

        """
        table_id = self.cursor.execute(
            f"SELECT {table}Id FROM {table} WHERE {column}='{value}'").fetchone()
        return table_id[0] if table_id else table_id
