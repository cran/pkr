LinAUC = function(x, y) # down="Linear"
{
# Author: Kyun-Seop Bae k@acr.kr
# AUC and AUMC of one segment by the linear trapezoidal rule.
# Delegates to NonCompart::LinAUC for a single maintained engine.
# INPUT
#    x: time
#    y: concentration
# RETURNS
#    named vector: AUC, AUMC
  return(NonCompart::LinAUC(x, y))
}
