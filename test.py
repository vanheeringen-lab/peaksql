import PeakSQL


db = PeakSQL.database.DataBase()
data = {'danRer11': ['48h', '8somites', '80%epiboly', 'dome', 'shield'],
        'mm10':     ['early_2cell', '2cell', '4cell', '8cell',
                     'Mm_E7', 'Mm_E8', 'Mm_E9', 'Mm_E10', 'Mm_E12', 'Mm_E14', 'Mm_E16', 'Mm_E18'],
        'Spur_5.0': ['28h', '18h', '24h', '30h', '39h', '50h', '60h', '70h'],
        'oryLat2':  ['st11', 'st13', 'st15', 'st19', 'st21', 'st24', 'st25', 'st28',
                     'st32', 'st36', 'st40']}

data = {'danRer11': ['48h', '8somites'],
        'mm10':     ['early_2cell', '2cell']}

# for assembly, stages in data.items():
#     # if assembly not in db.assemblies:
#     db.add_assembly(f'/vol/genomes/{assembly}/{assembly}.fa', assembly=assembly)

    # for stage in stages:
    #     db.add_data(f'/vol/atac-seq/macs2/{assembly}-{stage}_peaks.narrowPeak', assembly, stage)

import time
import numpy as np

now = time.time()
for i in range(10000):
    low, high = np.random.randint(50), np.random.randint(100, 200)
    db.cursor.execute(f"SELECT substr(sequence, {low}, {high}) FROM Chromosome where ChromosomeId=3")
    print(db.cursor.fetchone())
print(time.time() - now)
