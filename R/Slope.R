Slope = function(x, y)
{
# Author: Kyun-Seop Bae k@acr.kr
# Simple linear regression used for terminal-slope estimation.
# Delegates to NonCompart::Slope for a single maintained engine.
# INPUT
#    x: time
#    y: natural log of concentration
# RETURNS
#    named vector: R2, R2ADJ, LAMZNPT, LAMZ, b0, CORRXY, LAMZLL, LAMZUL, CLSTP
  return(NonCompart::Slope(x, y))
}
