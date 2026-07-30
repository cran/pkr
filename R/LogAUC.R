LogAUC = function(x, y) # down="Log" means Linear-Up Log-Down
{
# Author: Kyun-Seop Bae k@acr.kr
# AUC and AUMC by the linear-up/log-down trapezoidal rule.
# Delegates to NonCompart::LogAUC for a single maintained engine.
# INPUT
#    x: time
#    y: concentration
# RETURNS
#    named vector: AUC, AUMC
  return(NonCompart::LogAUC(x, y))
}
