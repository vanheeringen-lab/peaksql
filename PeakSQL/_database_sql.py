"""
All the SQL stuff
"""
import os

import pybedtools
import pyfaidx


def add_assembly(self, assembly_file, assembly=None, species=None):
    """

    """
    # set defaults if none provided
    assembly = assembly if assembly is not None else os.path.basename(assembly_file).split('.')[0]
    species = species if species is not None else assembly

    # Make sure the assembly hasn't been added yet
    assert assembly not in self.assemblies, f"Assembly '{assembly}' has already been added " \
                                            f"to the database!"

    self.cursor.execute(f"INSERT INTO Species (assembly, name, abspath) "
                        f"VALUES ('{assembly}', '{species}', '{assembly_file}')")
    species_id = self.get_id('Species', 'Assembly', assembly)

    # now load the chromosome sizes
    fasta = pyfaidx.Fasta(assembly_file)
    self.open_files[assembly_file] = fasta
    for sequence_name in fasta.keys():
        sequence = fasta[sequence_name]
        self.cursor.execute(f"INSERT INTO Chromosome (SpeciesId, Size, Chromosome) "
                            f"    VALUES('{species_id}', '{len(sequence)}', '{sequence_name}')")

    self.conn.commit()


def add_data(self, data_file, assembly, condition=None):
    """

    """
    # check for supported filetype
    *_, extension = os.path.splitext(data_file)
    assert extension in ['.narrowPeak'], f"The file extension you choose is not supported, TODO"

    # check if species it belongs to has already been added to the database
    species_id = self.get_id('Species', 'Assembly', assembly)
    assert species_id, f"Assembly '{assembly}' has not been added to the database yet."

    # check if the condition does not already exist for this species
    self.cursor.execute(f"SELECT Condition FROM Condition "
                        f"WHERE  SpeciesId='{species_id}' "
                        f"AND    Condition='{condition}'").fetchone()
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
