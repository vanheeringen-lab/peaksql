import sqlite3
import os
from functools import lru_cache

import pyfaidx

import peaksql.tables as tables


class DataBase:
    """
    The DataBase class serves as an easy interface to store and retrieve NGS data in a
    peaksql database.
    """

    def __init__(self, db: str = "PeakSQL.sqlite", in_memory: bool = False):
        """
        :param db: the name (path) of the database, will create a new database if it
            doesn't exist yet.
        :param in_memory: whether to load a (pre-existing) database into memory. This
            can be useful in combination with one of the dataloaders for
            faster queries.
        """
        self.db = db

        # connect, and set a relatively high timeout number for multiprocessing
        self.conn = sqlite3.connect(db, timeout=30)
        self.cursor = self.conn.cursor()

        self.in_memory = in_memory
        if in_memory:
            # start a connection with our memory and move our database there
            dest = sqlite3.connect(":memory:")
            self.conn.backup(dest)

            # replace the old connection and cursor with our new in-memory connection
            self.conn = dest
            self.cursor = self.conn.cursor()

        # register all the tables (Assembly, Chromosome, Condition, Peak)
        for table in [table for table in dir(tables) if not table.startswith("__")]:
            virtual = "VIRTUAL" if "virtual" in getattr(tables, table).lower() else ""
            self.cursor.execute(
                f"CREATE {virtual} TABLE IF NOT EXISTS {getattr(tables, table)}"
            )

        self.conn.commit()

        # if we are loading a pre-existing database connect to all the assemblies
        self.cursor.execute("SELECT Assembly, AbsPath FROM Assembly")
        self.fastas = {
            assembly: pyfaidx.Fasta(abspath)
            for assembly, abspath in self.cursor.fetchall()
        }

    @lru_cache()
    def get_assembly_id(self, assembly_name: str) -> int:
        """
        Get the AssemblyId based on Assembly (name).

        :param assembly_name: name of the assembly
        :return: id of the assembly
        """
        result = self.cursor.execute(
            f"SELECT AssemblyId FROM Assembly "
            f"    WHERE Assembly='{assembly_name}' "
            f"LIMIT 1"
        ).fetchone()
        if result:
            return result[0]
        raise ValueError(f"Assembly {assembly_name} is not present in the database")

    @lru_cache(maxsize=2 ** 16)
    def get_chrom_id(self, assembly_id: int, chrom_name: str) -> int:
        """
        Get the ChromosomeId based on assemblyId and Chromosome (name).

        :param assembly_id: id of the assembly
        :param chrom_name: name of the assembly
        :return: id of the chromosome
        """
        result = self.cursor.execute(
            f"SELECT ChromosomeId FROM Chromosome "
            f"    WHERE Chromosome='{chrom_name}' "
            f"    AND AssemblyId='{assembly_id}'"
            f"LIMIT 1"
        ).fetchone()
        if result:
            return result[0]
        raise ValueError(
            f"No chromosome {chrom_name} for assembly with assembly id {assembly_id}"
        )

    @property
    def assemblies(self):
        """
        All assemblies registred in the database.
        """
        return [
            val[0]
            for val in self.cursor.execute(f"SELECT Assembly FROM Assembly").fetchall()
        ]

    def add_assembly(
        self, assembly_path: str, assembly: str = None, species: str = None
    ):
        """
        Add an assembly (genome) to the database. Sequences from the assembly are
        retrieved with PyFaidx, so they aren't stored in the database, only the path to
        the fasta file is stored. This thus assumes that the assembly does not change
        location during during the database's lifetime.

        :param assembly_path: The path to the assembly file.
        :param assembly: The name of the assembly (optional: default is the name of the
            file).
        :param species: The name of the species the assembly belongs to (optional:
            default is the assembly name)
        """
        assert not self.in_memory, (
            "It is currently not supported to add data with an in-memory " "database."
        )
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

        fasta = pyfaidx.Fasta(abs_path)
        size = sum(len(seq) for seq in fasta.values())
        # add the assembly to the assembly table
        self.cursor.execute(
            f"INSERT INTO Assembly (Assembly, Species, Abspath, Size) "
            f"VALUES ('{assembly}', '{species}', '{abs_path}', {size})"
        )
        assembly_id = self.cursor.lastrowid

        # now fill the chromosome table
        offset = self.cursor.execute(
            f"SELECT SUM(Size) FROM Assembly WHERE AssemblyId < {assembly_id}"
        ).fetchone()[0]
        offset = 0 if offset is None else offset
        for sequence_name, sequence in fasta.items():
            size = len(sequence)
            self.cursor.execute(
                f"INSERT INTO Chromosome (AssemblyId, Size, Chromosome, Offset) "
                f"    VALUES('{assembly_id}', '{size}', '{sequence_name}', '{offset}')"
            )
            offset += size

        # clean up after yourself
        self.conn.commit()

    def add_data(self, data_path: str, assembly: str, condition: str = None):
        """
        Add data (bed or narrowPeak) to the database. The files are stored line by line

        :param data_path: The path to the assembly file.
        :param assembly: The name of the assembly. Requires the assembly to be added to
            the database prior.
        :param condition: Experimental condition (optional). This allows for filtering
            on conditions , e.g. when streaming data with a DataSet.
        """
        assert not self.in_memory, (
            "It is currently not supported to add data with an in-memory " "database."
        )
        # check for supported filetype
        *_, extension = os.path.splitext(data_path)
        # TODO: add more extensions
        assert extension in [
            ".narrowPeak",
            ".bed",
            ".bw",
        ], f"The file extension you choose is not supported"

        # check if species it belongs to has already been added to the database
        assembly_id = self.cursor.execute(
            f"SELECT AssemblyId FROM Assembly WHERE Assembly='{assembly}'"
        ).fetchone()
        assembly_id = assembly_id[0] if assembly_id else assembly_id
        assert assembly_id, (
            f"Assembly '{assembly}' has not been added to the database yet. Before "
            f"adding data you should add assemblies with the DataBase.add_assembly "
            f"method."
        )

        # Make sure that condition 'None' exists
        # This somehow locks the database when in __init__
        self.cursor.execute(
            "INSERT INTO Condition(ConditionId, Condition) SELECT 0, NULL "
            "WHERE NOT EXISTS(SELECT * FROM Condition WHERE ConditionId = 0)"
        )

        # get the condition id
        condition_id = (
            self.cursor.execute(
                f"SELECT ConditionId FROM Condition WHERE Condition='{condition}'"
            ).fetchone()
            if condition
            else (0,)
        )
        condition_id = condition_id[0] if condition_id else None

        # add the condition if necessary
        if condition and not condition_id:
            self.cursor.execute(f"INSERT INTO Condition VALUES(NULL, '{condition}')")

        if extension in [".bed", ".narrowPeak"]:
            self._add_bed(data_path, assembly_id, condition_id, extension)
        elif extension in [".bw"]:
            self._add_bigwig(data_path, assembly_id, condition_id, extension)

    def _add_bed(self, data_path, assembly_id, condition_id, extension):
        bed_lines = []
        virt_lines = []

        # get the current BedId we are at
        highest_id = self.cursor.execute(
            "SELECT BedId FROM Bed ORDER BY BedId DESC LIMIT 1"
        ).fetchone()
        if not highest_id:
            highest_id = 0
        else:
            highest_id = highest_id[0]

        with open(data_path) as bedfile:
            for i, line in enumerate(bedfile):
                bed = line.strip().split("\t")
                # print(bed)
                chromosome_id = self.get_chrom_id(assembly_id, bed[0])
                offset = self.cursor.execute(
                    f"SELECT Offset FROM Chromosome "
                    f"WHERE ChromosomeId = {chromosome_id}"
                ).fetchone()[0]

                chromstart, chromend = bed[1:3]
                chromstart = int(chromstart) + offset
                chromend = int(chromend) + offset
                virt_lines.append((i + 1 + highest_id, chromstart, chromend))

                if len(bed) == 3 and extension == ".bed":
                    bed_lines.append((condition_id, chromosome_id, None))
                elif len(bed) == 10 and extension == ".narrowPeak":
                    bed_lines.append((condition_id, chromosome_id, bed[9]))
                else:
                    fields = {".bed": 3, ".narrowPeak": 10}
                    assert False, (
                        f"Extension {extension} should have {fields[extension]} fields,"
                        f" however it has {len(fields)}"
                    )

        self.cursor.executemany(f"INSERT INTO Bed VALUES(NULL, ?, ?, ?)", bed_lines)

        # also add each bed entry to the BedVirtual table
        self.cursor.executemany(
            f"INSERT INTO BedVirtual VALUES(?, ?, ?)", virt_lines,
        )

        self.conn.commit()

    def _add_bigwig(self, data_path, assembly_id, condition_id, extension):
        raise NotImplementedError
