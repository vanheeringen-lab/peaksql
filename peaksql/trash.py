# @staticmethod
# def get_locations(chromstart, chromend):
#     return slice(chromstart, chromend)
#
# def get_peak_locations(self) -> dict:
#     """
#     Map of all assemblies and genomic coordinates and whether or not it falls within a peak
#     region.
#
#     :return
#     {
#     assembly1:
#         {chrom1: np.array([False, False, True, ..., True], dtype=bool),
#          chrom2: np.array([False, True, False, ..., True], dtype=bool),
#          ...
#     }
#     """
#     # raise NotImplementedError
#     total = 0
#     peaks = dict()
#     # loop over all the fetched results, where chromendsummit is either chromend or the summit
#     process = self._get_process()
#     for assembly, chrom, chromstart, chromendsummit in self.fetchall:
#         if assembly not in peaks:
#             peaks[assembly] = dict()
#         if chrom not in peaks[assembly]:
#             peaks[assembly][chrom] = np.zeros(
#                 len(self.databases[process].fastas[assembly][chrom]),
#                 dtype=bool
#             )
#             total += peaks[assembly][chrom].nbytes
#         # we get either a slice
#         peaks[assembly][chrom][self.get_locations(chromstart, chromendsummit)] = True
#
#     print("TOTAL", total, process)
#     return peaks
#
# def estimate_label_occurence(self, ) -> (int, int):
#     """
#     Get the number of occurrences of each label.
#
#     First returns the number of negatives (non-peak), and then the number of positives (peak).
#     """
#     raise NotImplementedError
#     peaks = 0
#     for i, (assembly, chrom) in enumerate(self.chromosomes[1:]):
#         for chromstart in self.positions[i + 1]:
#             peaks += at_least(self.peaks[assembly][chrom][chromstart:chromstart+self.seq_length]
#                               , self.frac_region)
#
#     return len(self) - peaks, peaks
