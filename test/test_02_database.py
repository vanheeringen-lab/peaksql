import unittest
import sys

import peaksql


DATABASE_BED = "test_peaksql_bed.sqlite"
DATABASE_NWP = "test_peaksql_narrowpeak.sqlite"


class TestDataBase(unittest.TestCase):
    """ A test class to test the peaksql.DATABASE_BED """

    def setUp(self):
        """
        Due to Unit Test architecture it will give us lots of ResourceWarnings which
        aren't relevant to our tests. We turn them off.
        """
        if not sys.warnoptions:
            import warnings

            warnings.simplefilter("ignore")

    def test_201_load_genomes(self):
        db = peaksql.DataBase(DATABASE_BED)
        db.add_assembly("test/data/assembly1.fa")
        db.add_assembly("test/data/assembly2.fa")
        self.assertEquals(db.assemblies, ["assembly1", "assembly2"])
        self.assertRaises(AssertionError, db.add_assembly, "test/data/assembly1.fa")

    def test_202_test_get_assembly_id(self):
        db = peaksql.DataBase(DATABASE_BED)
        self.assertEquals(db.get_assembly_id("assembly1"), 1)
        self.assertEquals(db.get_assembly_id("assembly2"), 2)
        self.assertEquals(db.get_assembly_id("assembly2"), 2)
        self.assertRaises(ValueError, db.get_assembly_id, "assembly3")

    def test_203_get_chrom_id(self):
        db = peaksql.DataBase(DATABASE_BED)
        self.assertEquals(db.get_chrom_id(1, "chr1"), 1)
        self.assertEquals(db.get_chrom_id(1, "chr2"), 2)
        self.assertEquals(db.get_chrom_id(2, "chr1"), 3)
        self.assertEquals(db.get_chrom_id(2, "chr3"), 4)
        self.assertRaises(ValueError, db.get_chrom_id, 1, 3)

    def test_204_add_bed(self):
        db = peaksql.DataBase(DATABASE_BED)
        db.add_data("test/data/assembly1.bed", "assembly1")
        db.add_data("test/data/assembly2.bed", "assembly2")
        assert db.cursor.execute("SELECT COUNT(BedId) FROM BED").fetchone() == (3,)
        assert db.cursor.execute(
            "SELECT Chromosome, ChromStart, ChromEnd "
            "FROM BED B "
            "JOIN Chromosome C on B.ChromosomeId = C.ChromosomeId "
            "JOIN Assembly A on C.AssemblyId = A.AssemblyId "
            "WHERE Assembly='assembly2' and Chromosome='chr1'"
        ).fetchall() == [("chr1", 20, 40)]

    def test_205_add_narrowpeak(self):
        db = peaksql.DataBase(DATABASE_NWP)
        db.add_assembly("test/data/assembly1.fa")
        db.add_data("test/data/assembly1.narrowPeak", "assembly1")

    def test_206_in_memory(self):
        db_file = peaksql.DataBase(DATABASE_BED, in_memory=False)
        db_memo = peaksql.DataBase(DATABASE_BED, in_memory=True)

        assert (
            db_file.cursor.execute("SELECT * FROM Assembly").fetchall()
            == db_memo.cursor.execute("SELECT * FROM Assembly").fetchall()
        )
        assert (
            db_file.cursor.execute("SELECT * FROM Chromosome").fetchall()
            == db_memo.cursor.execute("SELECT * FROM Chromosome").fetchall()
        )
        assert (
            db_file.cursor.execute("SELECT * FROM BED").fetchall()
            == db_memo.cursor.execute("SELECT * FROM BED").fetchall()
        )
