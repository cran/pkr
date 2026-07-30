Unit = function(code="", timeUnit="h", concUnit="ng/mL", doseUnit="mg", MW=0)
{
# Author: Kyun-Seop Bae k@acr.kr
# Unit strings and internal conversion factors for NCA parameters.
# Delegates to NonCompart::Unit for a single maintained engine.
# INPUT
#    code: SDTM PPTESTCD ("" returns all codes)
#    timeUnit: time unit
#    concUnit: concentration unit
#    doseUnit: dose unit (not amount per kg or per m2)
#    MW: molecular weight
# RETURNS
#    data.frame (or one row) of Unit and Factor, row-named by PPTESTCD
  return(NonCompart::Unit(code=code, timeUnit=timeUnit, concUnit=concUnit, doseUnit=doseUnit, MW=MW))
}
