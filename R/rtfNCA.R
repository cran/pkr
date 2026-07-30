rtfNCA = function(fileName="Temp-NCA.rtf", concData, colSubj="Subject", colTime="Time", colConc="conc", dose=0, adm="Extravascular", dur=0, doseUnit="mg", timeUnit="h", concUnit="ug/L", down="Linear", MW=0)
{
  if (!requireNamespace("rtf", quietly=TRUE)) stop("Package 'rtf' is needed for rtfNCA(). Please install it with install.packages('rtf').")

  nDev0 = length(dev.list())
  rtf = rtf::RTF(fileName)
# Ensure the RTF file is finalized, and any device rtf::addPlot() opened for a
# failed panel is closed, even if a subject errors mid-loop.
  on.exit({ try(rtf::done(rtf), silent=TRUE); while (length(dev.list()) > nDev0) dev.off() })
  rtf::addHeader(rtf, title="Individual Noncompartmental Analysis Result")
  rtf::addNewLine(rtf)
  rtf::addHeader(rtf, "Table of Contents")
  rtf::addTOC(rtf)
  rtf::setFontSize(rtf, font.size=10)

  tAll = as.numeric(concData[,colTime])
  cAll = as.numeric(concData[,colConc])
  tAll = tAll[is.finite(tAll)]
  cAll = cAll[is.finite(cAll)]
  if (length(tAll) == 0 | length(cAll) == 0) stop("No finite time or concentration value!")
  maxx = max(tAll)
  maxy = max(cAll)
  cPos = cAll[cAll > 0]
  miny = if (length(cPos) > 0) min(cPos) else maxy

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

    rtf::addPageBreak(rtf)
    rtf::addHeader(rtf, paste("Subject ID =", cID), TOC.level=1)
    for (j in seq_along(tRes)) rtf::addParagraph(rtf, tRes[j])

    rtf::addPageBreak(rtf)
    rtf::addHeader(rtf, paste("Subject ID =", cID))
    rtf::addPlot(rtf, plot.fun=plot, width=6, height=4, res=300, x=x, y=y, type="b", cex=0.7,
          xlim=c(0,maxx), ylim=c(0,maxy),
          xlab=paste0("Time (", timeUnit, ")"), ylab=paste0("Concentration (", concUnit, ")"))
    rtf::addPlot(rtf, plot.fun=Plot4rtf, width=6, height=4, res=300, x=x, y=y, type="b", cex=0.7,
          xlim=c(0, maxx), ylim=c(miny, maxy),
          xlab=paste0("Time (", timeUnit, ")"), ylab=paste0("Concentration (log interval) (", concUnit, ")"), tabRes=tabRes)
  }

  rtf::done(rtf)
}
