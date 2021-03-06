readEX = function(folders)
{
  EX = combXPT(folders, "EX")

# Check DoseUnit
  DoseUnit =  unique(Trim(EX[,"EXDOSU"]))
  DoseUnit = DoseUnit[nchar(DoseUnit) > 0]
  if (length(DoseUnit) != 1) warning("Dose unit is missing or various!")
   
# Dosing unit should not be composite one like mg/kg, mg/m2
  if (length(strsplit(DoseUnit, "/")[[1]]) != 1) stop("Dose unit should not be based on body weight or BSA!")     

# If EXENDTC is empty, set it as EXSTDTC.
  for (i in 1:nrow(EX)) {
    if (EX[i,"EXENDTC"] == "") EX[i,"EXENDTC"] = EX[i,"EXSTDTC"]
  }

# Numeric type column will be set
  colNum = intersect(c("EXDOSE", "EXTPTNUM"), colnames(EX))
  nCol = length(colNum)
  for (i in 1:nCol) {
    EX[,colNum[i]] = as.double(EX[,colNum[i]])
  }

  if (!("EXROUTE" %in% colnames(EX))) EX[,"EXROUTE"] = "ORAL" # assume oral administration

  return(EX)
}
