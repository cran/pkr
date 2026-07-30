readEX = function(folders)
{
  EX = combXPT(folders, "EX")

# Check DoseUnit
  DoseUnit =  unique(Trim(EX[,"EXDOSU"]))
  DoseUnit = DoseUnit[nchar(DoseUnit) > 0]
  if (length(DoseUnit) != 1) warning("Dose unit is missing or various!")

# Dosing unit should not be composite one like mg/kg, mg/m2
# (guard so a missing/various unit warns above instead of erroring here)
  if (any(lengths(strsplit(DoseUnit, "/")) != 1)) stop("Dose unit should not be based on body weight or BSA!")

# If EXENDTC is empty, set it as EXSTDTC.
  if (!("EXENDTC" %in% colnames(EX))) EX[,"EXENDTC"] = ""
  emptyEnd = is.na(EX[,"EXENDTC"]) | Trim(EX[,"EXENDTC"]) == ""
  EX[emptyEnd, "EXENDTC"] = EX[emptyEnd, "EXSTDTC"]

# Numeric type column will be set
  colNum = intersect(c("EXDOSE", "EXTPTNUM"), colnames(EX))
  for (i in seq_along(colNum)) {
    EX[,colNum[i]] = as.double(EX[,colNum[i]])
  }

  if (!("EXROUTE" %in% colnames(EX))) EX[,"EXROUTE"] = "ORAL" # assume oral administration

  return(EX)
}
