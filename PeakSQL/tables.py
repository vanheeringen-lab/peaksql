"""
Collection of all the tables used by PeakSQL.
"""
# Species table
spe = "Assembly (" \
      "    AssemblyId INTEGER PRIMARY KEY AUTOINCREMENT," \
      "    Assembly TEXT NOT NULL," \
      "    Species," \
      "    AbsPath TEXT UNIQUE NOT NULL" \
      ")"

# Chromosome table
chr = "Chromosome (" \
      "    ChromosomeId INTEGER PRIMARY KEY AUTOINCREMENT," \
      "    Chromosome TEXT," \
      "    Size INT NOT NULL," \
      "    Sequence TEXT NOT NULL," \
      "    AssemblyId NOT NULL," \
      "    FOREIGN KEY(AssemblyId) REFERENCES Assembly(AssemblyId)" \
      ")"

# Condition table
con = "Condition (" \
      "    ConditionId INTEGER PRIMARY KEY AUTOINCREMENT, " \
      "    Condition TEXT," \
      "    AssemblyId NOT NULL," \
      "    FOREIGN KEY(AssemblyId) REFERENCES Assebly(AssemblyId)" \
      ")"

# Bed table
bed = ("Bed (" 
       "    BedId INTEGER PRIMARY KEY AUTOINCREMENT, " 
       "    AssemblyId NOT NULL," 
       "    ConditionId NOT NULL," 
       "    ChromosomeId NOT NULL,"    # 1.  chrom (chrom ID not name) (REQUIRED)
       "    ChromStart INT NOT NULL,"  # 2.  chromStart                (REQUIRED)
       "    ChromStop INT NOT NULL,"   # 3.  chromEnd                  (REQUIRED)
       "    Name TEXT,"                # 4.  name
       "    Score INT,"                # 5.  score
       "    Strand TEXT,"              # 6.  strand
       "    SignalValue TEXT,"         # 7.  signalValue
       "    PValue FLOAT,"             # 8.  pValue
       "    QValue FLOAT,"             # 9.  qValue
       "    Peak INT,"                 # 10. peak (chromStart + Peak)
       "    FOREIGN KEY(AssemblyId)   REFERENCES Assembly(AssemblyId)," 
       "    FOREIGN KEY(ChromosomeId) REFERENCES Chromosome(ChromosomeId)," 
       "    FOREIGN KEY(ConditionId)  REFERENCES Condition(ConditionId)" 
       ")")
