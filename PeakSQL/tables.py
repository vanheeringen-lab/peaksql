"""
Collection of all the tables used by PeakSQL.
"""
# Species table
spe = "Species (" \
      "    SpeciesId INTEGER PRIMARY KEY AUTOINCREMENT," \
      "    Assembly TEXT UNIQUE," \
      "    Name TEXT" \
      ")"

# Chromosome table
chr = "Chromosome (" \
      "    ChromosomeId INTEGER PRIMARY KEY AUTOINCREMENT," \
      "    Chromosome TEXT," \
      "    Size INT NOT NULL," \
      "    SpeciesId," \
      "    FOREIGN KEY(SpeciesId) REFERENCES Species(SpeciesId)" \
      ")"

# Condition table
con = "Condition (" \
      "    ConditionId INTEGER PRIMARY KEY AUTOINCREMENT, " \
      "    Condition TEXT," \
      "    SpeciesId," \
      "    FOREIGN KEY(SpeciesId) REFERENCES Species(SpeciesId)" \
      ")"

# Bed table
bed = ("Bed (" 
       "    BedId INTEGER PRIMARY KEY AUTOINCREMENT, " 
       "    SpeciesId," 
       "    ConditionId," 
       "    ChromosomeId,"             # 1.  chrom (chrom ID not name) (REQUIRED)
       "    ChromStart INT NOT NULL,"  # 2.  chromStart                (REQUIRED)
       "    ChromStop INT NOT NULL,"   # 3.  chromEnd                  (REQUIRED)
       "    Name TEXT,"                # 4.  name
       "    Score INT,"                # 5.  score
       "    Strand TEXT,"              # 6.  strand
       "    SignalValue TEXT,"         # 7.  signalValue
       "    PValue FLOAT,"             # 8.  pValue
       "    QValue FLOAT,"             # 9.  qValue
       "    Peak INT,"                 # 10. peak (chromStart + Peak)
       "    FOREIGN KEY(SpeciesId)    REFERENCES Species(SpeciesId)," 
       "    FOREIGN KEY(ChromosomeId) REFERENCES Chromosome(ChromosomeId)," 
       "    FOREIGN KEY(ConditionId)  REFERENCES Condition(ConditionId)" 
       ")")
