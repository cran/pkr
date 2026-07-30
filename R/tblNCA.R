tblNCA = function(concData, key="Subject", colTime="Time", colConc="conc", dose=0,
         adm="Extravascular", dur=0, doseUnit="mg", timeUnit="h", concUnit="ug/L",
         down="Linear", MW=0, returnNA=FALSE)
{
# Author: Kyun-Seop Bae k@acr.kr
# Do NCA for every subject (key) and return a result table.
# Delegates to NonCompart::tblNCA so that pkr and the author's NonCompart package
# share one maintained engine. NonCompart::tblNCA already handles multiple keys
# cleanly (the key columns keep their own names) and binds rows by name, so the
# earlier positional-rbind and multi-key naming problems no longer arise.
# INPUT
#   concData: concentration data table
#   key: column name(s) identifying each profile
#   colTime: column name for time
#   colConc: column name for concentration
#   dose: dose (scalar or one per profile)
#   adm: "Extravascular", "Bolus", or "Infusion"
#   dur: duration of infusion
#   doseUnit, timeUnit, concUnit: units
#   down: "Linear" or "Log"
#   MW: molecular weight
#   returnNA: if FALSE, drop columns that are entirely NA across all profiles
# RETURNS
#   data.frame of NCA results, with a "units" attribute

  Res = NonCompart::tblNCA(concData, key=key, colTime=colTime, colConc=colConc,
                           dose=dose, adm=adm, dur=dur, doseUnit=doseUnit,
                           timeUnit=timeUnit, concUnit=concUnit, down=down, MW=MW)

  if (!isTRUE(returnNA)) {
    Units = attr(Res, "units")
    keep = vapply(seq_len(ncol(Res)), function(j) !all(is.na(Res[[j]])), logical(1))
    keep[seq_along(key)] = TRUE  # always keep the key columns
    Res2 = Res[, keep, drop=FALSE]
    if (!is.null(Units)) attr(Res2, "units") = Units[keep]
    Res = Res2
  }
  return(Res)
}
