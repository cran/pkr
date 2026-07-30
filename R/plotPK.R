plotPK = function(concData, id, Time, conc, unitTime="hr", unitConc="ng/mL",
                  trt="", fit="Linear", dose=0, adm="Extravascular", dur=0,
                  outdir="Output", name="")
{
# Author: Kyun-Seop Bae k@acr.kr
# Write the standard set of PK figures for a concentration data table:
#   individual linear and semi-log profiles (one PDF each) and pooled linear,
#   semi-log and mean +- 95% CI profiles (one raster image each).
# INPUT
#   concData: concentration data table
#   id, Time, conc: column names for subject, time and concentration
#   unitTime, unitConc: axis units
#   trt: column name for treatment, or "" for a single-treatment study
#   fit: "Linear" or "Log" for the downward trapezoidal rule
#   dose, adm, dur: dosing information passed on to NCA()
#   outdir: directory for the figures; "" means the working directory
#   name: stem used in the file names; "" means the deparsed concData argument
# RETURNS
#   invisibly, the paths of the files written

  if (!is.numeric(dose) | !is.numeric(dur) | !is.character(adm) | !is.character(fit)) stop("Bad Input!")
  for (cName in c(id, Time, conc, if (nzchar(trt)) trt)) {
    if (!(cName %in% colnames(concData))) stop(paste0("There is no column named '", cName, "'!"))
  }

  SUBJIDs = unique(as.character(concData[,id]))
  nSUBJID = length(SUBJIDs)
  if (nSUBJID == 0) stop("No subject in concData!")

  if (nzchar(trt)) {
    TRTs = sort(unique(as.character(concData[,trt])))
  } else {
    TRTs = ""
  }
  nTRT = length(TRTs)

  nCell = nSUBJID*nTRT
  if (length(dose) > 1 & length(dose) != nCell) stop("dose should be fixed or given for each subject!")
  if (length(dur) > 1 & length(dur) != nCell) stop("dur should be fixed or given for each subject!")

  if (!nzchar(name)) name = .pkSafeName(deparse(substitute(concData)))
  DrugName = name

# Axis ranges. Non-finite values are dropped rather than propagated: a single NA
# concentration used to make every ylim non-finite and abort the whole function.
  tv = as.numeric(concData[,Time])
  cv = as.numeric(concData[,conc])
  tv = tv[is.finite(tv)]
  cv = cv[is.finite(cv)]
  if (length(tv) == 0 | length(cv) == 0) stop("No finite time or concentration value!")
  timemin = min(tv) ; timemax = max(tv)
  concmin = min(cv) ; concmax = max(cv)
  cvPos = cv[cv > 0]
  if (length(cvPos) == 0) {
    warning("No positive concentration; the semi-log figures are skipped.")
    logAx = NULL
  } else {
    logAx = .pkLogAxis(min(cvPos), max(cvPos))
  }

  Result = NCA(concData, id, Time, conc, trt=trt, fit=fit, dose=dose, adm=adm,
               dur=dur, uTime=unitTime, uConc=unitConc)

  if (!nzchar(outdir)) outdir = "."
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(outdir)) stop(paste0("Could not create the output directory '", outdir, "'!"))

# Every device opened below is closed again even if a panel errors; otherwise a
# half-written file stays current and the caller's next plot lands inside it.
  nDev0 = length(dev.list())
  on.exit(while (length(dev.list()) > nDev0) dev.off())

  subCols = .pkPalette(nSUBJID)
  subPchs = .pkPch(nSUBJID)
  trtCols = .pkPalette(nTRT)
  written = character(0)

  cellDat = function(i, j) {
    keep = as.character(concData[,id]) == SUBJIDs[i]
    if (nzchar(trt)) keep = keep & as.character(concData[,trt]) == TRTs[j]
    concData[which(keep), , drop=FALSE]
  }
  cellRes = function(i, j) {
    keep = as.character(Result[,id]) == SUBJIDs[i]
    if (nzchar(trt)) keep = keep & as.character(Result[,trt]) == TRTs[j]
    Result[which(keep), , drop=FALSE]
  }
  panelTitle = function(i, j) {
    if (nzchar(trt)) paste("Subject ID", SUBJIDs[i], TRTs[j]) else paste("Subject ID", SUBJIDs[i])
  }
  noteNCA = function(Res, left, right, col) {
    val = function(nm) if (nm %in% colnames(Res) && nrow(Res) > 0) signif(as.numeric(Res[1,nm]), 3) else NA
    mtext(paste0(left[1], ": ", val(left[2]), " ", left[3]), side=3, adj=0, col=col, cex=0.9)
    mtext(paste0(right[1], ": ", val(right[2]), " ", right[3]), side=3, adj=1, col=col, cex=0.9)
  }
  blankPanel = function(msg="") {
    plot.new()
    if (nzchar(msg)) text(0.5, 0.5, msg, cex=0.9, col="grey40")
  }

# ---------------------------------------------------------------- individual
# One panel per subject and treatment. Panels are laid out in a grid capped at
# four columns; a single row of nTRT panels used to make both the page width and
# the margin request grow without bound.
  individual = function(logScale)
  {
    if (logScale & is.null(logAx)) return(invisible(NULL))
    label = if (logScale) "Individual PK Log 10 Scale for " else "Individual PK Linear Scale for "
    fName = file.path(outdir, paste0(label, DrugName, ".pdf"))
    if (nzchar(trt)) {
      g = .pkGrid(nTRT)
    } else {
      g = c(2, 2)
    }
    nr = g[1] ; nc = g[2]
    pdf(file=fName, height=min(5*nr, 20), width=min(5*nc, 20))
    on.exit(dev.off())
    par(oma=c(0, 0, 2, 0), mfrow=c(nr, nc), mar=c(4.5, 4.5, 3.5, 1.5))
    for (i in seq_len(nSUBJID)) {
      for (j in seq_len(nTRT)) {
        Dat = cellDat(i, j)
        x = as.numeric(Dat[,Time])
        y = as.numeric(Dat[,conc])
        ok = if (logScale) is.finite(x) & is.finite(y) & y > 0 else is.finite(x) & is.finite(y)
        cCol = if (nzchar(trt)) trtCols[j] else subCols[1]
        if (sum(ok) == 0) {
          blankPanel(paste0(panelTitle(i, j), "\n(no plottable point)"))
          next
        }
        x = x[ok] ; y = y[ok]
        if (logScale) {
          plot(x, log10(y), pch=16, type="b", col=cCol, yaxt="n",
               xlab=paste0("Time (", unitTime, ")"),
               ylab=paste0("Concentration (", unitConc, ")"),
               ylim=logAx$ylim, xlim=c(timemin, timemax))
          .pkDrawLogAxis(logAx)
        } else {
          plot(x, y, pch=16, type="b", col=cCol,
               xlab=paste0("Time (", unitTime, ")"),
               ylab=paste0("Concentration (", unitConc, ")"),
               ylim=c(concmin, concmax), xlim=c(timemin, timemax))
        }
        title(panelTitle(i, j), cex.main=1.1)
        Res = cellRes(i, j)
        if (logScale) {
          noteNCA(Res, c("Tmax", "TMAX", unitTime), c("Half-life", "LAMZHL", unitTime), cCol)
        } else {
          noteNCA(Res, c("Cmax", "CMAX", unitConc),
                  c("AUClast", "AUCLST", paste(unitTime, unitConc, sep="*")), cCol)
        }
      }
# Pad each subject out to a full page so that pages stay aligned with subjects
# even when a treatment cell is empty.
      if (nzchar(trt) & nTRT < nr*nc) for (k in seq_len(nr*nc - nTRT)) blankPanel()
    }
    mtext(paste("Individual Profiles of", DrugName), outer=TRUE, cex=1.2)
    written <<- c(written, fName)
    return(invisible(NULL))
  }

# -------------------------------------------------------------------- pooled
# All subjects overlaid, one panel per treatment, with the subject legend in its
# own layout cell so that it can never overlap the curves or spill off the page.
  pooled = function(logScale)
  {
    if (logScale & is.null(logAx)) return(invisible(NULL))
    label = if (logScale) "PK Profile Log 10 Scale for " else "PK Profile Linear Scale for "
    fName = file.path(outdir, paste0(label, DrugName, ".tiff"))
    g = .pkGrid(nTRT)
    nr = g[1] ; nc = g[2]
    .pkRaster(fName, nr, nc, panel=5, legend=1.9, res=150)
    on.exit(dev.off())
    par(oma=c(0, 0, 3, 0))
    lay = cbind(matrix(seq_len(nr*nc), nrow=nr, ncol=nc, byrow=TRUE), nr*nc + 1)
    layout(lay, widths=c(rep(5, nc), 1.9))
    par(mar=c(4.5, 4.5, 3, 1.5))
    for (j in seq_len(nr*nc)) {
      if (j > nTRT) { blankPanel() ; next }
      if (logScale) {
        plot(NA, NA, yaxt="n", xlab=paste0("Time (", unitTime, ")"),
             ylab=paste0("Concentration (", unitConc, ")"),
             ylim=logAx$ylim, xlim=c(timemin, timemax))
        .pkDrawLogAxis(logAx)
      } else {
        plot(NA, NA, xlab=paste0("Time (", unitTime, ")"),
             ylab=paste0("Concentration (", unitConc, ")"),
             ylim=c(concmin, concmax), xlim=c(timemin, timemax))
      }
      if (nzchar(trt)) title(TRTs[j], col.main="indianred", cex.main=1.1)
      for (i in seq_len(nSUBJID)) {
        Dat = cellDat(i, j)
        x = as.numeric(Dat[,Time])
        y = as.numeric(Dat[,conc])
        ok = if (logScale) is.finite(x) & is.finite(y) & y > 0 else is.finite(x) & is.finite(y)
        if (sum(ok) == 0) next
        points(x[ok], if (logScale) log10(y[ok]) else y[ok], pch=subPchs[i], type="b",
               col=subCols[i], cex=0.8)
      }
    }
    .pkLegendPanel(SUBJIDs, col=subCols, pch=subPchs, title="Subject ID")
    mtext(paste("Concentration vs. Time Profile of", DrugName), outer=TRUE, cex=1.3)
    written <<- c(written, fName)
    return(invisible(NULL))
  }

# ------------------------------------------------------------- mean with CI
  meanCI = function()
  {
    fName = file.path(outdir, paste0("PK Profile with CI for ", DrugName, ".tiff"))
    g = .pkGrid(nTRT)
    nr = g[1] ; nc = g[2]
    .pkRaster(fName, nr, nc, panel=5, res=150)
    on.exit(dev.off())
    par(oma=c(0, 0, 3, 0), mfrow=c(nr, nc), mar=c(5.5, 4.5, 4, 1.5))
    anyPanel = FALSE
    for (j in seq_len(nTRT)) {
      DatT = if (nzchar(trt)) concData[which(as.character(concData[,trt]) == TRTs[j]), , drop=FALSE] else concData
      nom = .pkNominalTime(as.numeric(DatT[,Time]), nrow(cellDat(1, j)))
      nDropped = 0
      Dat3 = NULL
      for (i in seq_len(nSUBJID)) {
        Dat = cellDat(i, j)
        if (nrow(Dat) == 0) next
        if (nrow(Dat) != length(nom$t)) { nDropped = nDropped + 1 ; next }
        Dat$NomTime = nom$t
        Dat3 = rbind(Dat3, Dat)
      }
      if (nDropped > 0) {
        warning(paste0(nDropped, " subject(s) were left out of the mean profile",
                       if (nzchar(trt)) paste0(" of ", TRTs[j]) else "",
                       " because their sampling count differs from the nominal schedule."))
      }
      if (is.null(Dat3) || nrow(Dat3) == 0) {
        blankPanel(paste0(if (nzchar(trt)) paste0(TRTs[j], "\n") else "",
                          "no subject matches the\nnominal sampling schedule"))
        next
      }
      anyPanel = TRUE
      meanC = aggregate(as.numeric(Dat3[,conc]), by=list(Dat3$NomTime), FUN=function(v) mean(v, na.rm=TRUE))
      sdC   = aggregate(as.numeric(Dat3[,conc]), by=list(Dat3$NomTime), FUN=function(v) sd(v, na.rm=TRUE))
      nC    = aggregate(as.numeric(Dat3[,conc]), by=list(Dat3$NomTime), FUN=NROW)
      tt = meanC[,1]
      mm = meanC$x
      err = qnorm(0.975)*sdC$x/sqrt(nC$x)
      err[!is.finite(err)] = 0
      upper = mm + err
      lower = mm - err
# The error bars, not just the means, decide the y range; with the fixed
# c(concmin, concmax) the upper whiskers were clipped off the panel.
      yr = range(c(mm, upper, lower, concmin, concmax), na.rm=TRUE, finite=TRUE)
      halfCap = diff(range(c(tt, timemin, timemax)))*0.01
      plot(tt, mm, col=trtCols[j], pch=16,
           xlab=paste0("Time (", unitTime, ")"),
           ylab=paste0("Concentration (", unitConc, ")"),
           ylim=yr, xlim=c(timemin, timemax))
      if (nzchar(trt)) title(TRTs[j], col.main="olivedrab", cex.main=1.1)
      segments(tt, lower, tt, upper, col="grey")
      segments(tt - halfCap, lower, tt + halfCap, lower, col="grey")
      segments(tt - halfCap, upper, tt + halfCap, upper, col="grey")
      legend("topright", c("Mean", "95% CI"), pch=c(16, NA), lty=c(NA, 1),
             col=c(trtCols[j], "grey"), inset=0.05, cex=0.7)
# paste() is vectorised, so passing the whole count vector to mtext() drew every
# bin count on top of the others at the same line. The note sits in the bottom
# margin, where a long list of bins cannot run into the legend.
      mtext(paste("Counts in each bin:", paste(nom$ct, collapse=", ")), side=1, line=4, adj=0, cex=0.6)
      Res = if (nzchar(trt)) Result[which(as.character(Result[,trt]) == TRTs[j]), , drop=FALSE] else Result
      mCmax = if ("CMAX" %in% colnames(Res)) signif(mean(as.numeric(Res[,"CMAX"]), na.rm=TRUE), 3) else NA
      mAUC  = if ("AUCLST" %in% colnames(Res)) signif(mean(as.numeric(Res[,"AUCLST"]), na.rm=TRUE), 3) else NA
      mtext(paste0("Mean Cmax: ", mCmax, " ", unitConc), side=3, adj=0, col=trtCols[j], cex=0.7)
      mtext(paste0("Mean AUClast: ", mAUC, " ", unitTime, "*", unitConc), side=3, adj=1,
            col=trtCols[j], cex=0.7)
    }
    if (nTRT < nr*nc) for (k in seq_len(nr*nc - nTRT)) blankPanel()
    mtext(paste("Concentration vs. Time Profile of", DrugName), outer=TRUE, cex=1.3)
    if (anyPanel) written <<- c(written, fName)
    return(invisible(NULL))
  }

  individual(logScale=FALSE)
  individual(logScale=TRUE)
  pooled(logScale=FALSE)
  pooled(logScale=TRUE)
  meanCI()

  return(invisible(written))
}
