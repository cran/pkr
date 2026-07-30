AUC = function(x, y, down="Linear")
{
# Author: Kyun-Seop Bae k@acr.kr
# Cumulative AUC and AUMC. Delegates to NonCompart::AUC so that pkr and the
# author's NonCompart package share a single, maintained implementation.
# INPUT
#    x: time or similar vector
#    y: concentration or similar vector
#    down: "Linear" or "Log" for the downward trapezoidal rule
# RETURNS
#    n x 2 matrix of cumulative AUC and AUMC
  return(NonCompart::AUC(x, y, down=down))
}
