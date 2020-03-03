"""
Collection of all the tables used by PeakSQL.
"""
# Assembly table
ASS = (
    "Assembly ("
    "    AssemblyId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    Assembly TEXT NOT NULL,"
    "    Species,"
    "    Size INT NOT NULL,"
    "    AbsPath TEXT UNIQUE NOT NULL"
    ")"
)

# Chromosome table
CHR = (
    "Chromosome ("
    "    ChromosomeId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    Chromosome TEXT,"
    "    Size INT NOT NULL,"
    "    Offset INT NOT NULL,"
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
BED = (
    "Bed ("
    "    BedId INTEGER PRIMARY KEY AUTOINCREMENT,"
    "    ConditionId,"
    "    ChromosomeId NOT NULL,"
    "    Peak INT,"  # narrowPeak summit
    "    FOREIGN KEY(ChromosomeId) REFERENCES Chromosome(ChromosomeId),"
    "    FOREIGN KEY(ConditionId)  REFERENCES Condition(ConditionId)"
    ")"
)

# Virtual Bed table, complement of the BED table. Uses r*tree for faster queries
BED_VIRT = (
    f"BedVirtual USING rtree("
    f"    BedId INT, "  # Foreign key not implemented :(
    f"    ChromStart INT, "
    f"    ChromEnd INT, "
    f")"
)
