readPC = function(folders)
{
# Author: Kyun-Seop Bae k@acr.kr
# Read and combine PC domain files, keeping plasma, done records; BLQ -> 0.
  PC = combXPT(folders, "PC")

# Provide defaults for optional SDTM columns that the code below relies on.
  if (!("PCSPEC" %in% colnames(PC))) PC[,"PCSPEC"] = "PLASMA"
  if (!("PCLLOQ" %in% colnames(PC))) PC[,"PCLLOQ"] = "0"
  if (!("PCSTAT" %in% colnames(PC))) PC[,"PCSTAT"] = ""
  if (!("PCORRES" %in% colnames(PC))) PC[,"PCORRES"] = ""

# Treat below-limit-of-quantitation original results as 0 (before numeric coercion).
  codeBLQ = c("< 0", "<0", "NQ", "BLQ", "BQL", "BQOL", "BLOQ", "<LOQ", "<LLOQ")
  PC[UT(PC[,"PCORRES"]) %in% codeBLQ, "PCSTRESN"] = 0

# Numeric type column will be set
  colNum = intersect(c("PCSTRESN", "VISITNUM", "PCTPTNUM"), colnames(PC))
  for (i in seq_along(colNum)) {
    PC[,colNum[i]] = as.double(PC[,colNum[i]])
  }

  PC = PC[UT(PC[,"PCSPEC"]) == "PLASMA",] # currently PLASMA only
  PC = PC[UT(PC[,"PCSTAT"]) != "NOT DONE",] # remove not done

  return(PC)
}
