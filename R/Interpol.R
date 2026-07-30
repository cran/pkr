Interpol = function(x, y, xnew, Slope=0, b0=0, down="Linear")
{
# Author: Kyun-Seop Bae k@acr.kr
# Insert an interpolated (x, y) point into a concentration-time series.
# Delegates to NonCompart::Interpol for a single maintained engine.
# INPUT
#   x: time
#   y: concentration
#   xnew: new time point to interpolate
#   Slope: terminal slope of log(y) ~ x
#   b0: intercept of log(y) ~ x
#   down: "Linear" or "Log"
# RETURNS
#   list(x, y) with the new point inserted
  return(NonCompart::Interpol(x, y, xnew, Slope=Slope, b0=b0, down=down))
}
