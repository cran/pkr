sNCA = function(x, y, dose=0, adm="Extravascular", dur=0, doseUnit="mg", timeUnit="h", concUnit="ug/L", iAUC="", down="Linear", MW=0, returnNA=TRUE, UsePoints=NULL)
{
# Author: Kyun-Seop Bae k@acr.kr
# Single-subject noncompartmental analysis.
# The computation is delegated to NonCompart::sNCA so that pkr and the author's
# NonCompart package always produce identical results from one maintained engine.
# NonCompart::sNCA returns only the parameters relevant to the administration
# route (dropping the NA ones). When returnNA=TRUE (the historical pkr default)
# the result is expanded to the full, fixed 48-parameter vector, padding
# route-irrelevant parameters with NA, so downstream pkr code and existing user
# scripts see a stable column set and order.
# INPUT
#   x: time vector
#   y: concentration vector
#   dose: given dose (not per kg or per m2)
#   adm: "Extravascular", "Bolus", or "Infusion"
#   dur: duration of infusion
#   doseUnit, timeUnit, concUnit: units
#   iAUC: data.frame(Name, Start, End) for interval AUCs, or "" for none
#   down: "Linear" or "Log"
#   MW: molecular weight
#   returnNA: if TRUE, return the full 48-parameter vector padded with NA
#   UsePoints: indices into x and y, after non-finite points are dropped, to be
#     used for the terminal slope instead of the automatic search. NULL leaves the
#     automatic selection in place.
# RETURNS
#   named numeric vector of NCA parameters, with a "units" attribute

  Res = NonCompart::sNCA(x, y, dose=dose, adm=adm, dur=dur, doseUnit=doseUnit,
                         timeUnit=timeUnit, concUnit=concUnit, iAUC=iAUC,
                         down=down, MW=MW, UsePoints=UsePoints)

  if (!isTRUE(returnNA)) return(Res)
  UsedPoints = attr(Res, "UsedPoints")

# Expand to the full, fixed pkr parameter set (returnNA=TRUE contract).
  RetNames1 = c("b0", "CMAX", "CMAXD", "TMAX", "C0", "AUCPBEO", "AUCPBEP", "TLAG", "CLST", "CLSTP", "TLST", "LAMZHL", "LAMZ",
             "LAMZLL", "LAMZUL", "LAMZNPT", "CORRXY", "R2", "R2ADJ", "AUCLST", "AUCALL",
             "AUCIFO", "AUCIFOD", "AUCIFP", "AUCIFPD", "AUCPEO", "AUCPEP",
             "AUMCLST", "AUMCIFO", "AUMCIFP", "AUMCPEO", "AUMCPEP",
             "MRTIVLST", "MRTIVIFO", "MRTIVIFP", "MRTEVLST", "MRTEVIFO", "MRTEVIFP",
             "VZO", "VZP", "VZFO", "VZFP", "CLO", "CLP", "CLFO", "CLFP", "VSSO", "VSSP")
  iAUCNames = setdiff(names(Res), RetNames1) # any interval-AUC columns appended by NonCompart
  fullNames = c(RetNames1, iAUCNames)

  Out = rep(NA_real_, length(fullNames))
  names(Out) = fullNames
  Out[names(Res)] = as.numeric(Res)

  Units = NonCompart::Unit(doseUnit=doseUnit, timeUnit=timeUnit, concUnit=concUnit, MW=MW)
  attr(Out, "units") = c(Units[RetNames1, 1], rep(Units["AUCLST", 1], length(iAUCNames)))
# Which points the terminal slope came from is worth keeping through the padding.
  if (!is.null(UsedPoints)) attr(Out, "UsedPoints") = UsedPoints
  return(Out)
}
