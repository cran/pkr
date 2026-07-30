pdfNCA = function(fileName="Temp-NCA.pdf", concData, colSubj="Subject", colTime="Time", colConc="conc", dose=0, adm="Extravascular", dur=0, doseUnit="mg", timeUnit="h", concUnit="ug/L", down="Linear", MW=0) 
{
# par(no.readonly=TRUE) opens the default device when none is current, which used
# to leave a stray Rplots.pdf behind. Only the PDF opened here is touched, and it
# is closed again below, so nothing of the caller's needs saving or restoring.
  nDev0 = length(dev.list())
  PrepPDF(fileName)
# Ensure the PDF device is closed even if a subject errors mid-loop.
  on.exit(while (length(dev.list()) > nDev0) dev.off())
  AddPage()
  Text1(1, 1, "Individual Noncompartmental Analysis Result", Cex=1.2)

  tAll = as.numeric(concData[,colTime])
  cAll = as.numeric(concData[,colConc])
  tAll = tAll[is.finite(tAll)]
  cAll = cAll[is.finite(cAll)]
  if (length(tAll) == 0 | length(cAll) == 0) stop("No finite time or concentration value!")
  maxx = max(tAll)
  maxy = max(cAll)
  cPos = cAll[cAll > 0]
  miny = if (length(cPos) > 0) min(cPos) else NA
# One shared decade-aligned axis for every subject, so the semi-log panels can
# be compared page to page. The tick range used to be taken from each subject
# while ylim came from the whole study, which could leave the axis unlabelled.
  logAx = if (is.na(miny)) NULL else .pkLogAxis(miny, maxy)

  IDs = unique(concData[,colSubj])
  nID = length(IDs)
  Res = vector()
  for (i in seq_len(nID)) {
    cID = IDs[i]
    x = concData[concData[,colSubj]==cID, colTime]
    y = concData[concData[,colSubj]==cID, colConc]
    tabRes = sNCA(x, y, dose=dose, adm=adm, dur=dur, doseUnit=doseUnit, timeUnit=timeUnit, concUnit=concUnit, down=down, MW=MW)
    tRes = txtNCA(x, y, dose=dose, adm=adm, dur=dur, doseUnit=doseUnit, timeUnit=timeUnit, concUnit=concUnit, down=down, MW=MW, returnNA=FALSE)
    Res = c(Res, tRes)

    AddPage(Header1=paste("Subject ID =", cID))
    TextM(tRes, StartRow=1, Header1=paste("Subject ID =", cID))

    scrnmat = matrix(0, 3, 4)
    scrnmat[1,] = c(0, 1, 0, 1)
    scrnmat[2,] = c(0.1, 0.9, 0.50, 0.95)
    scrnmat[3,] = c(0.1, 0.9, 0.05, 0.50)
    ScrNo = split.screen(scrnmat)

    screen(ScrNo[1])
    par(adj=0)
    Text1(1, 1, paste("Subject ID =", cID), Cex=1.0)

    screen(ScrNo[2])
    par(oma=c(1,1,1,1), mar=c(4,4,3,1), adj=0.5)
    plot(x, y, type="b", cex=0.7, xlim=c(0,maxx), ylim=c(0,maxy), xlab=paste0("Time (", timeUnit, ")"), ylab=paste0("Concentration (", concUnit, ")"))

    screen(ScrNo[3])
    par(oma=c(1,1,1,1), mar=c(4,4,3,1), adj=0.5)
    ok = is.finite(x) & is.finite(y) & y > 0
    x0 = x[ok]
    y0 = y[ok]
    if (length(y0) > 0 & !is.null(logAx)) {  # a subject with no positive concentration has no semi-log panel
      plot(x0, log10(y0), type="b", cex=0.7, xlim=c(0, maxx), ylim=logAx$ylim, yaxt="n", xlab=paste0("Time (", timeUnit, ")"), ylab=paste0("Concentration (log interval) (", concUnit, ")"))
      .pkDrawLogAxis(logAx)
      if (isTRUE(is.finite(tabRes["LAMZ"]))) {
        x1 = tabRes["LAMZLL"]
        x2 = tabRes["LAMZUL"]
        deltaX = x1 * 0.05
        y1 = log10(exp(1))*(tabRes["b0"] - tabRes["LAMZ"] * (x1 - deltaX))
        y2 = log10(exp(1))*(tabRes["b0"] - tabRes["LAMZ"] * (x2 + deltaX))
        lines(c(x1 - deltaX, x2 + deltaX), c(y1, y2), lty=2, col="red")
      }
    }
    close.screen(all.screens=TRUE)
  }
# device close and par restore handled by on.exit()
}
