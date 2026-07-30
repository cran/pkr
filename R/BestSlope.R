BestSlope = function(x, y, adm="Extravascular", TOL=1e-4)
{
# Author: Kyun-Seop Bae k@acr.kr
# Automatic terminal-slope (lambda z) selection by the maximum adjusted R-squared
# criterion. Delegates to NonCompart::BestSlope for a single maintained engine.
# INPUT
#    x: time
#    y: concentration (original scale, not logged)
#    adm: administration method, "Extravascular", "Bolus", or "Infusion"
#    TOL: tolerance for adjusted R-squared comparison
# RETURNS
#    named vector: R2, R2ADJ, LAMZNPT, LAMZ, b0, CORRXY, LAMZLL, LAMZUL, CLSTP
  return(NonCompart::BestSlope(x, y, adm=adm, TOL=TOL))
}
