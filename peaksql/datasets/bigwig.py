# import numpy as np
import pyBigWig

from .basedataset import _DataSet


class BigWigDataSet(_DataSet):
    """
    The BigWigDataSet ...
    """

    def __init__(
        self, database: str, where: str = "", seq_length: int = 200, **kwargs: int
    ):
        _DataSet.__init__(self, database, where, seq_length, **kwargs)
        self.open = pyBigWig.open(
            "/vol/joost/macs2/rheMac10-Moorlag-H3K27ac-NHP-LvE-1-17867.bw"
        )

    # @profile
    def get_label(self, assembly, chrom, chromstart, chromend):
        # assemblyid = self.database.get_assembly_id(assembly)
        # chromosomeid = self.database.get_chrom_id(assemblyid, chrom)

        # print(chrom, chromstart, chromend)
        # try:
        labels = self.open.values(chrom, chromstart, chromend)
        # except:
        #     labels = 0
        # bed_virtual = f"BedVirtual_{assemblyid}"
        # query = f"""
        #     SELECT {self.SELECT_LABEL} FROM {bed_virtual}
        #     INNER JOIN Bed on {bed_virtual}.BedId = Bed.BedId
        #     WHERE ({chromstart} <= {bed_virtual}.ChromEnd) AND
        #           ({chromend} >= {bed_virtual}.ChromStart)
        # """
        # query_result = self.database.cursor.execute(query).fetchall()
        #
        # positions = self.array_from_query(
        #     query_result, chromosomeid, chromstart, chromend
        # )
        # labels = self.label_from_array(positions)

        return labels
