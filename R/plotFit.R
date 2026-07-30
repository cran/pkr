plotFit = function(concData, id, Time, conc, mol="", adm="Extravascular", ID="", Mol="")
{
# Author: Kyun-Seop Bae k@acr.kr
# Plot one subject's concentration-time profile on a natural-log axis with the
# selected terminal slope drawn over it.
# INPUT
#   concData: concentration data table
#   id, Time, conc: column names for subject, time and concentration
#   mol: column name for the analyte, or "" if there is only one
#   adm: "Extravascular", "Bolus", or "Infusion"
#   ID: subject to plot
#   Mol: analyte to plot, used only when mol is given
# RETURNS
#   the BestSlope() result for that subject

  keep = as.character(concData[,id]) == as.character(ID)
  if (nzchar(mol)) keep = keep & as.character(concData[,mol]) == as.character(Mol)
  keep[is.na(keep)] = FALSE
  x = as.numeric(concData[which(keep), Time])
  y = as.numeric(concData[which(keep), conc])
  if (length(x) == 0) stop(paste0("No record for ID '", ID, "'", if (nzchar(mol)) paste0(" and ", Mol) else "", "!"))

  finalMat = BestSlope(x, y, adm=adm)

# Only positive concentrations exist on a log axis. The old code substituted 0.1
# for zeros, which put an arbitrary point on the plot and into the axis range.
  ok = is.finite(x) & is.finite(y) & y > 0
  if (sum(ok) == 0) stop("No positive concentration to plot!")
  xp = x[ok]
  yp = log(y[ok])

# The fitted intercept is included in the range so that the slope line is never
# drawn off the panel, and so that an extreme b0 cannot reverse the limits.
  yr = range(c(yp, if (is.finite(finalMat["b0"])) finalMat["b0"]), finite=TRUE)
  yr = c(floor(yr[1]), ceiling(yr[2]))
  if (yr[2] <= yr[1]) yr[2] = yr[1] + 1

  plot(xp, yp, yaxt="n", ylim=yr, xlab="Time", ylab="Concentration",
       main=paste("Best Fit ID:", ID))
  yticks = seq(yr[1], yr[2], by=max(1, ceiling(diff(yr)/9)))
  axis(2, at=yticks, labels=as.expression(lapply(yticks, function(i) bquote(e^.(i)))))
  if (isTRUE(finalMat["LAMZNPT"] > 0) && is.finite(finalMat["b0"]) && is.finite(finalMat["LAMZ"])) {
    abline(a=finalMat["b0"], b=-finalMat["LAMZ"], untf=TRUE, col="blue")
  }

  return(finalMat)
}
