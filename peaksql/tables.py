"""
Collection of all the tables used by PeakSQL.
"""
# Assembly table
ASS = (
    "Assembly ("
    "    AssemblyId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    Assembly TEXT NOT NULL,"
    "    Species,"
    "    AbsPath TEXT UNIQUE NOT NULL"
    ")"
)

# Chromosome table
CHR = (
    "Chromosome ("
    "    ChromosomeId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    Chromosome TEXT,"
    "    Size INT NOT NULL,"
    "    AssemblyId NOT NULL,"
    "    FOREIGN KEY(AssemblyId) REFERENCES Assembly(AssemblyId)"
    ")"
)

# Condition table
CON = (
    "Condition ("
    "    ConditionId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    Condition TEXT"
    ")"
)

# Bed table
# fmt: off
BED = (
    "Bed ("
    "    BedId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    ConditionId,"
    "    ChromosomeId NOT NULL,"    # 1.  chrom (chrom ID not name) (REQUIRED)
    "    ChromStart INT NOT NULL,"  # 2.  chromStart                (REQUIRED)
    "    ChromEnd INT NOT NULL,"    # 3.  chromEnd                  (REQUIRED)
    "    Name TEXT,"                # 4.  name
    "    Score INT,"                # 5.  score
    "    Strand TEXT,"              # 6.  strand
    "    SignalValue TEXT,"         # 7.  signalValue
    "    PValue FLOAT,"             # 8.  pValue
    "    QValue FLOAT,"             # 9.  qValue
    "    Peak INT,"                 # 10. peak (chromStart + Peak)
    "    FOREIGN KEY(ChromosomeId) REFERENCES Chromosome(ChromosomeId),"
    "    FOREIGN KEY(ConditionId)  REFERENCES Condition(ConditionId)"
    ")"
)
# fmt: on

# BigWig table
WIG = (
    "Wig ("
    "    WigID INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    AbsPath TEXT UNIQUE NOT NULL,"
    "    ConditionId INT NOT NULL,"
    "    AssemblyId INT NOT NULL,"
    "    FOREIGN KEY(ConditionId) REFERENCES Condition(ConditionId),"
    "    FOREIGN KEY(AssemblyId) REFERENCES Assembly(AssemblyId)"
    ")"
)
