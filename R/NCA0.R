NCA0 = function(EX0, PC0, fit="Linear", excludeInfusion=TRUE, label="")
{
# Author: Kyun-Seop Bae k@acr.kr
# Noncompartmental analysis of one dosing record and its concentration records.
# INPUT
#   EX0: one EX record. A named character vector, a one-row data.frame or a
#        one-row matrix are all accepted.
#   PC0: that record's PC rows, with PCDTC, PCSTRESN and PCSTRESU
#   fit: "Linear" or "Log" for the downward trapezoidal rule
#   excludeInfusion: for an infusion, estimate the terminal slope from
#     post-infusion samples only; the default. Has no effect on a bolus or an
#     extravascular dose, which have no infusion to exclude.
#     The restriction lives here rather than in sNCA() because this is the layer
#     that knows the dosing record, and so the infusion end time. sNCA() stays a
#     faithful wrapper around NonCompart::sNCA.
#   label: identifier used only in warning messages
# RETURNS
#   named character vector: AdmMethod, Dose, UnitDose, Dur, UnitTime, UnitConc,
#   then the NCA parameters

# Every access below is by name, so a one-row data.frame - which is what
# subsetting a data.frame by row actually yields, and what rNCA passed for years -
# gave a one-column data.frame instead of a string and failed downstream inside
# NonCompart::Unit(). Accept the shape and flatten it here.
  if (!is.null(dim(EX0))) EX0 = unlist(lapply(as.data.frame(EX0, stringsAsFactors=FALSE)[1,], as.character))

  cDose = as.numeric(EX0["EXDOSE"])
  cUnitDose = as.character(EX0["EXDOSU"])
  cRoute = toupper(Trim(as.character(EX0["EXROUTE"])))
  cStart = as.character(EX0["EXSTDTC"])
  cEnd = as.character(EX0["EXENDTC"])
# Route decides the model; the recorded duration only distinguishes an infusion
# from a bolus WITHIN the intravenous route. Testing the duration first, as this
# used to, classified any extravascular record carrying an end date/time - routine
# in real EX domains - as an infusion, and so reported it with VZO/CLO/VSS in
# place of VZFO/CLFO/MRTEV.
  hasEnd = !is.na(cEnd) & cEnd != "NA" & nzchar(cEnd) & !is.na(cStart) & cEnd > cStart
  if (cRoute == "INTRAVENOUS") {
    if (hasEnd) {
      admMethod = "Infusion"
      cDur = as.numeric(difftime(strptime(cEnd,"%Y-%m-%dT%H:%M:%S"), strptime(cStart, "%Y-%m-%dT%H:%M:%S"), units="hours"))
    } else {
      admMethod = "Bolus"
      cDur = 0
    }
  } else {
    admMethod = "Extravascular"
    cDur = 0
  }
  cUnitConc = unique(Trim(as.character(PC0[,"PCSTRESU"])))
  cUnitConc = cUnitConc[nchar(cUnitConc) > 0]
# sNCA takes one concentration unit; passing a vector silently used only the
# first while the unit table was built from all of them.
  if (length(cUnitConc) > 1) {
    warning(paste0("Various concentration units (", paste(cUnitConc, collapse=", "), "); the first is used."))
  }
  cUnitConc = if (length(cUnitConc) == 0) "" else cUnitConc[1]
  TAD = as.numeric(difftime(strptime(PC0[,"PCDTC"],"%Y-%m-%dT%H:%M:%S"), strptime(cStart, "%Y-%m-%dT%H:%M:%S"), units="hours"))
  TAD[TAD < 0] = 0
  y = as.numeric(PC0[,"PCSTRESN"])

  UsePoints = NULL
  if (isTRUE(excludeInfusion) && admMethod == "Infusion") {
    UsePoints = .postInfusionPoints(TAD, y, cDur, admMethod)
    if (is.null(UsePoints)) {
# Falling back silently would report a terminal slope drawn from the very samples
# the caller asked to exclude, so say what happened.
      warning(paste0("Too few post-infusion samples to estimate the terminal slope",
                     if (nzchar(label)) paste0(" for ", label) else "",
                     "; the automatic selection over all samples is used instead."))
    }
  }

  Res0 = sNCA(TAD, y, dose=cDose, adm=admMethod, dur=cDur, doseUnit=cUnitDose, concUnit=cUnitConc, down=fit, returnNA=FALSE, UsePoints=UsePoints)
  tNames = names(Res0)
  Res = c(admMethod, cDose, cUnitDose, cDur, "h", cUnitConc, Res0)
  names(Res) = c("AdmMethod", "Dose", "UnitDose", "Dur", "UnitTime", "UnitConc", tNames)
  return(Res)
}
