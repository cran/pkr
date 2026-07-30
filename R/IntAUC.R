IntAUC = function(x, y, t1, t2, Res, down="Linear")
{
# Author: Kyun-Seop Bae k@acr.kr
# Interval (partial) AUC between t1 and t2, interpolating end points as needed.
# Delegates to NonCompart::IntAUC for a single maintained engine.
# INPUT
#    x: time
#    y: concentration
#   t1: start time (interpolated if not in x)
#   t2: end time (interpolated if not in x)
#  Res: result of sNCA (provides TLST, LAMZ, b0)
#  down: "Linear" or "Log"
# RETURNS
#    interval AUC (numeric)
  return(NonCompart::IntAUC(x, y, t1, t2, Res, down=down))
}
