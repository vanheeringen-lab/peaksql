import PeakSQL


db = PeakSQL.DataBase()
for assembly in ['mm10', 'danRer11', 'dm6']:
    if assembly not in db.assemblies:
        db.add_assembly(f'/vol/genomes/{assembly}/{assembly}.fa')

db.add_data('/vol/atac-seq/macs2/mm10-DRR138928_peaks.narrowPeak', 'mm10', 'test')
