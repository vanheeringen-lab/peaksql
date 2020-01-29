"""

"""
import sqlite3
import os

import pybedtools
import pyfaidx

import peaksql.tables as tables


class DataBase:
    """

    """

    def __init__(self, db="PeakSQL.sqlite"):
        self.db = db

        self.conn = sqlite3.connect(db)
        self.cursor = self.conn.cursor()

        self.open_files = {}

        # register all the tables (Assembly, Chromosome, Condition, Peak)
        for table in [table for table in dir(tables) if not table.startswith("__")]:
            self.cursor.execute(f"CREATE TABLE IF NOT EXISTS {getattr(tables, table)}")

        self.conn.commit()

        # if we are loading pre-existing database
        self.cursor.execute("SELECT Assembly, AbsPath FROM Assembly")
        self.fastas = {
            assembly: pyfaidx.Fasta(abspath)
            for assembly, abspath in self.cursor.fetchall()
        }

    @property
    def assemblies(self):
        """
        return all registred assemblies
        """
        return [
            val[0]
            for val in self.cursor.execute(f"SELECT Assembly FROM Assembly").fetchall()
        ]

    def add_assembly(self, assembly_file, assembly=None, species=None):
        """

        """
        # set defaults if none provided
        assembly = (
            assembly if assembly else os.path.basename(assembly_file).split(".")[0]
        )
        species = species if species else assembly

        # TODO: what should default behaviour be (probably overwrite)?
        if assembly in self.assemblies:
            return

        # Make sure the assembly hasn't been added yet
        assert (
            assembly not in self.assemblies
        ), f"Assembly '{assembly}' has already been added to the database!"

        self.cursor.execute(
            f"INSERT INTO Assembly (Assembly, Species, Abspath) "
            f"VALUES ('{assembly}', '{species}', '{assembly_file}')"
        )

        # TODO: walrus operator in python3.8
        assembly_id = self.cursor.execute(
            f"SELECT AssemblyId FROM Assembly WHERE Assembly='{assembly}'"
        ).fetchone()[0]

        # now load the chromosome sizes
        fasta = pyfaidx.Fasta(assembly_file)
        self.open_files[assembly_file] = fasta
        for sequence_name, sequence in fasta.items():
            self.cursor.execute(
                f"INSERT INTO Chromosome (AssemblyId, Size, Chromosome) "
                f"    VALUES('{assembly_id}', '{len(sequence)}', '{sequence_name}')"
            )

        self.conn.commit()

    def add_data(self, data_file, assembly, condition=None):
        """

        """
        # check for supported filetype
        *_, extension = os.path.splitext(data_file)
        assert extension in [
            ".narrowPeak"
        ], f"The file extension you choose is not supported, TODO"

        # check if species it belongs to has already been added to the database
        # TODO: walrus operator in python3.8
        assembly_id = self.cursor.execute(
            f"SELECT AssemblyId FROM Assembly WHERE Assembly='{assembly}'"
        ).fetchone()
        assembly_id = assembly_id[0] if assembly_id else assembly_id
        assert (
            assembly_id
        ), f"Assembly '{assembly}' has not been added to the database yet."

        # check if the condition does not already exist for this species
        self.cursor.execute(
            f"SELECT Condition FROM Condition "
            f"WHERE  AssemblyId='{assembly_id}' "
            f"AND    Condition='{condition}'"
        ).fetchone()
        # TODO: what should default behaviour be (probably overwrite)?
        if self.cursor.fetchone() is not None:
            return
        assert not self.cursor.fetchone(), (
            f"Condition '{condition}' has already been added to "
            f"the database for assembly '{assembly}'."
        )

        # add the condition
        self.cursor.execute(
            f"INSERT INTO Condition VALUES(NULL, '{condition}', '{assembly_id}')"
        )

        # TODO: walrus operator in python3.8
        condition_id = self.cursor.execute(
            f"SELECT ConditionId FROM Condition WHERE Condition='{condition}'"
        ).fetchone()
        condition_id = condition_id[0] if condition_id else condition_id
        bed = pybedtools.BedTool(data_file)
        for region in bed:
            # TODO: walrus operator in python3.8
            chromosome_id = self.cursor.execute(
                f"SELECT ChromosomeId FROM Chromosome "
                f"    WHERE Chromosome='{region.chrom}' "
                f"    AND AssemblyId='{assembly_id}'"
            ).fetchone()
            chromosome_id = chromosome_id[0] if chromosome_id else chromosome_id

            # chromosome_id = self.get_id('Chromosome', 'Chromosome', region.chrom)
            if len(region.fields) == 3:
                raise NotImplementedError
            else:
                assert len(region.fields) == 10, "TODO error message"
                region.fields[-1] = int(region.fields[1]) + int(region.fields[-1])
                self.cursor.execute(
                    f"""INSERT INTO Bed """
                    f"""  VALUES(NULL, '{assembly_id}', '{condition_id}', """
                    f"""  '{chromosome_id}', '{"', '".join(region.fields[1:])}')"""
                )

        self.conn.commit()
