"""

"""
import sqlite3
import os
from functools import lru_cache

import pybedtools
import pyfaidx

import peaksql.tables as tables


class DataBase:
    """
    The DataBase class is an easy interface to store pre-processed NGS data for Machine Learning.

    It allows for ...
    """
    def __init__(self, db: str = "PeakSQL.sqlite", in_memory: bool = False):
        self.db = db

        # connect, and set a relatively high timeout number for multiprocessing
        self.conn = sqlite3.connect(db, timeout=30)
        self.cursor = self.conn.cursor()
        # WAL not really necessary here
        self.cursor.execute('PRAGMA journal_mode=WAL')

        self.in_memory = in_memory
        if in_memory:
            # start a connection with our memory and move our database there
            dest = sqlite3.connect(':memory:')
            self.conn.backup(dest)

            # replace the old connection and cursor with our new in-memory connection
            self.conn = dest
            self.cursor = self.conn.cursor()

        # register all the tables (Assembly, Chromosome, Condition, Peak)
        for table in [table for table in dir(tables) if not table.startswith("__")]:
            self.cursor.execute(f"CREATE TABLE IF NOT EXISTS {getattr(tables, table)}")

        self.conn.commit()

        # if we are loading a pre-existing database connect to all the assemblies
        self.cursor.execute("SELECT Assembly, AbsPath FROM Assembly")
        self.fastas = {
            assembly: pyfaidx.Fasta(abspath)
            for assembly, abspath in self.cursor.fetchall()
        }

    @lru_cache()
    def get_assembly_id(self, assembly_name):
        """
        Quickly get the AssemblyId based on Assembly (name).
        """
        return self.cursor.execute(
            f"SELECT AssemblyId FROM Assembly "
            f"    WHERE Assembly='{assembly_name}' "
            f"LIMIT 1"
        ).fetchone()[0]

    @lru_cache(maxsize=2**16)
    def get_chrom_id(self, assembly_id, chrom_name):
        """
        Quickly get the ChromosomeId based on assemblyId and Chromosome (name).
        """
        return self.cursor.execute(
            f"SELECT ChromosomeId FROM Chromosome "
            f"    WHERE Chromosome='{chrom_name}' "
            f"    AND AssemblyId='{assembly_id}'"
            f"LIMIT 1"
        ).fetchone()[0]

    @property
    def assemblies(self):
        """
        return all registred assemblies
        """
        return [
            val[0]
            for val in self.cursor.execute(f"SELECT Assembly FROM Assembly").fetchall()
        ]

    def add_assembly(self, assembly_path: str, assembly: str = None, species: str = None):
        """
        Add an assembly (genome) to the database. Sequences from the assembly are retrieved with
        PyFaidx, and the database assumes the assembly does not change location during during its
        lifetime.

        :param assembly_path: The path to the assembly file.
        :param assembly: The name of the assembly (optional: default is the name of the file).
        :param species: The name of the species the assembly belongs to (optional: default is the
        assembly name)
        """
        assert not self.in_memory, "It is currently not supported to add data with an in-memory " \
                                   "database."
        # set defaults if none provided
        assembly = (
            assembly if assembly else os.path.basename(assembly_path).split(".")[0]
        )
        species = species if species else assembly
        abs_path = os.path.abspath(assembly_path)

        # TODO: what should default behaviour be (probably overwrite) or ignore??
        # TODO: check for assembly_path instead of name?
        # make sure the assembly hasn't been added yet
        assert (
            assembly not in self.assemblies
        ), f"Assembly '{assembly}' has already been added to the database!"

        # add the assembly to the assembly table
        self.cursor.execute(
            f"INSERT INTO Assembly (Assembly, Species, Abspath) "
            f"VALUES ('{assembly}', '{species}', '{abs_path}')"
        )
        assembly_id = self.cursor.lastrowid

        # also add a virtual Bed table for R*Tree
        self.cursor.execute(f"CREATE VIRTUAL TABLE BedVirtual_{assembly_id} USING rtree("
                            f"    BedId INT, "  # Foreign key not implemented
                            f"    ChromStart INT, "
                            f"    ChromEnd INT, "
                            f")"
                            )

        # now fill the chromosome table
        fasta = pyfaidx.Fasta(abs_path)
        for sequence_name, sequence in fasta.items():
            self.cursor.execute(
                f"INSERT INTO Chromosome (AssemblyId, Size, Chromosome) "
                f"    VALUES('{assembly_id}', '{len(sequence)}', '{sequence_name}')"
            )

        # clean up after yourself
        self.conn.commit()

    def add_data(self, data_path: str, assembly: str, condition: str = None):
        """
        Add data (bed or narrowPeak) to the database. The files are stored line by line

        :param data_path: The path to the assembly file.
        :param assembly: The name of the assembly. Requires the assembly to be added to the database
        prior.
        :param condition: Experimental condition (optional). This allows for filtering on conditions
        , e.g. when streaming data with a DataSet.
        """
        assert not self.in_memory, "It is currently not supported to add data with an in-memory " \
                                   "database."
        # check for supported filetype
        *_, extension = os.path.splitext(data_path)
        # TODO: add more extensions
        assert extension in [
            ".narrowPeak",
            ".bed"
        ], f"The file extension you choose is not supported"

        # check if species it belongs to has already been added to the database
        assembly_id = self.cursor.execute(
            f"SELECT AssemblyId FROM Assembly WHERE Assembly='{assembly}'"
        ).fetchone()
        assembly_id = assembly_id[0] if assembly_id else assembly_id
        assert (
            assembly_id
        ), f"Assembly '{assembly}' has not been added to the database yet. Before adding data you" \
           f" should add assemblies with the DataBase.add_assembly method."

        # Make sure that condition 'None' exists
        self.cursor.execute("INSERT INTO Condition(ConditionId, Condition) SELECT 0, NULL "
                            "WHERE NOT EXISTS(SELECT * FROM Condition WHERE ConditionId = 0)")

        # get the condition id
        condition_id = self.cursor.execute(
            f"SELECT ConditionId FROM Condition WHERE Condition='{condition}'"
        ).fetchone()
        condition_id = condition_id[0] if condition_id else None

        # add the condition if necessary
        if condition and not condition_id:
            self.cursor.execute(
                f"INSERT INTO Condition VALUES(NULL, '{condition}')"
            )

        bed = pybedtools.BedTool(data_path)
        lines = []

        # get the current BedId we are at
        highest_id = self.cursor.execute("SELECT BedId FROM Bed ORDER BY BedId DESC LIMIT 1").fetchone()
        if not highest_id:
            highest_id = 0
        else:
            highest_id = highest_id[0]

        for region in bed:
            chromosome_id = self.get_chrom_id(assembly_id, region.chrom)

            if len(region.fields) == 3 and extension == '.bed':
                lines.append((condition_id, chromosome_id, *region.fields[1:3], None, None,
                              None, None, None, None, None))
            elif len(region.fields) == 10 and extension == '.narrowPeak':
                lines.append((condition_id, chromosome_id, *region.fields[1:]))
            else:
                fields = {'.bed': 3, '.narrowPeak': 10}
                assert False, f"Extension {extension} should have {fields[extension]} fields, " \
                              f"however it has {len(fields)}"

        self.cursor.executemany(
            f"""INSERT INTO Bed """
            f"""  VALUES(NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", lines
        )
        self.conn.commit()

        # TODO: reformat data loading
        # TODO: chromend and chromstart are twice in the database, once in Bed and once in
        # BedVirtual
        # also add each bed entry to the BedVirtual table
        for i, line in enumerate(lines):
            sth = f"INSERT INTO BedVirtual_{assembly_id}(BedId, ChromStart, ChromEnd) "\
                  f"VALUES({highest_id + i + 1}, {line[2]}, {line[3]})"
            self.cursor.execute(sth)

        self.conn.commit()

    def create_index(self, overwrite=False):
        """
        Create indexes for faster queries. It is highly recommended to run this before loading data.
        """
        if overwrite:
            # remove the index if exists and if we want to overwrite
            for index in ['Bed_ChromosomeId_Chromstart_Chromend_idx',
                        'Assembly_Assembly',
                        'Chromosome_Chromosome_AssemblyId']:
                self.cursor.execute(f"DROP INDEX IF EXISTS {index}")

        # exit when still an index exists, since then we do not want to overwrite
        indices = self.cursor.execute("PRAGMA INDEX_LIST('BED');").fetchall()
        indices_flat = [item for sublist in indices for item in sublist]
        if "Bed_ChromosomeId_Chromstart_Chromend_idx" in indices_flat:
            return

        # now add the indexes
        self.cursor.execute("CREATE INDEX Bed_ChromosomeId_Chromstart_Chromend_idx ON "
                            "BED (ChromosomeId, ChromStart, ChromEnd)")
        self.cursor.execute("CREATE INDEX Assembly_Assembly ON "
                            "Assembly (Assembly)")
        self.cursor.execute("CREATE INDEX Chromosome_Chromosome_AssemblyId ON "
                            "Chromosome (Chromosome, AssemblyId)")

        self.conn.commit()
