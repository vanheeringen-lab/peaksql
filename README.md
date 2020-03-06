# PeakSQL

[![Maintainability](https://api.codeclimate.com/v1/badges/d5f1443a164eb0d64d33/maintainability)](https://codeclimate.com/github/vanheeringen-lab/peaksql/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/d5f1443a164eb0d64d33/test_coverage)](https://codeclimate.com/github/vanheeringen-lab/peaksql/test_coverage)
![ContinuousIntegration](https://github.com/vanheeringen-lab/peaksql/workflows/ContinuousIntegration/badge.svg)
![ContinuousDeployment](https://github.com/vanheeringen-lab/peaksql/workflows/ContinuousDeployment/badge.svg)

Easy, fast, and dynamic machine learning database for genomics. Supports common bed-like dataformats like *.bed*, *.narrowPeak*, and *bedgraph*; and the binary *bigwig* format. 

### Installation
PeakSQL can be installed through pip:
```
pip install peaksql
```
Or installed from source:
```
git clone https://github.com/vanheeringen-lab/peaksql
cd peaksql
pip install .
```

### Getting started
```
import peaksql

# paths to our files
db_file = 'peakSQL.sqlite'  # where to store our database
assembly = "/path/to/hg19.fa"
data = "binding_sites.bed"

# load data into database
db = peaksql.database.DataBase(db_file)
db.add_assembly(assembly, assembly="hg19", species="human")
db.add_data(intervals_file, assembly="hg19")

# now load as dataset
dataset = peaksql.BedRegionDataSet(db_file, seq_length=101, stride=200)

# use the dataset in your application
for seq, label in dataset:
    ...
```
