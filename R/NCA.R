NCA = function(concData, id, Time, conc, trt="", fit="Linear", dose=0, adm="Extravascular", dur=0, report="Table", iAUC="", uTime="h", uConc="ug/L", uDose="mg")
{
# Author: Kyun-Seop Bae k@acr.kr
# Noncompartmental analysis for many subjects. Calls sNCA (delegated to
# NonCompart) for each subject/treatment.
  if (!is.numeric(dose) | !is.numeric(dur) | !is.character(adm) | !is.character(fit)) stop("Bad Input!")

  colOrd = paste0(adm, "Default")
  ColName00 = RptCfg[RptCfg[,colOrd] > 0, c("PPTESTCD", colOrd)]
  ColName0 = ColName00[order(ColName00[, colOrd]), "PPTESTCD"]

  if (!(max(dose) > 0)) ColName0 = setdiff(ColName0, c("CMAXD", "AUCIFOD", "AUCIFPD"))

  if (!missing(iAUC) && is.data.frame(iAUC)) {
    ColName0 = union(ColName0, as.character(iAUC[,"Name"]))
  }

  SUBJIDs = unique(as.character(concData[,id]))
  nSUBJID = length(SUBJIDs)

  if (trt == "") {
    TRTs = ""
    nTRT = 1
  } else {
    TRTs = sort(unique(as.character(concData[,trt])))
    nTRT = length(TRTs)
  }

  if (length(dose) > 1 & length(dose) != nSUBJID*nTRT) stop("dose should be fixed or given for each subject!")
  if (length(dur) > 1 & length(dur) != nSUBJID*nTRT) stop("dur should be fixed or given for each subject!")

# Pick the dose/dur for a given cell index (1-based over the SUBJID x TRT grid).
  pick = function(v, k) if (length(v) > 1) v[k] else v

  if (trt == "") {
    Res0 = data.frame(SUBJID=character(), stringsAsFactors=FALSE)
    if (report == "Table") Result = data.frame() else Result = vector()
    for (i in seq_len(nSUBJID)) {
      cSUBJID = SUBJIDs[i]
      Dat = concData[concData[,id]==cSUBJID,]
      if (nrow(Dat) > 0) {
        x = as.numeric(Dat[,Time])
        y = as.numeric(Dat[,conc])
        cDose = pick(dose, i)
        cTimeInfusion = pick(dur, i)
        if (adm == "Infusion" & !(cTimeInfusion > 0)) stop("Infusion mode should have dur larger than 0!")

        Res0 = rbind(Res0, data.frame(cSUBJID, stringsAsFactors=FALSE))
        cRes = sNCA(x, y, dose=cDose, adm=adm, dur=cTimeInfusion, iAUC=iAUC, doseUnit=uDose, timeUnit=uTime, concUnit=uConc, down=fit, returnNA=TRUE)
        cResult = cRes[ColName0]  # select by name in report order; fixed length per subject
        if (report == "Table") {
          Result = rbind(Result, cResult)
        } else {
          Result = c(Result, "NCA REPORT", paste0("Subject=", cSUBJID), "", cResult, "", "")
        }
      }
    }
    if (report == "Table") {
      Result = cbind(Res0, Result)
      colnames(Result) = c(id, ColName0)
      rownames(Result) = NULL
    }
  } else {
    Res0 = data.frame(SUBJID=character(), TRT=character(), stringsAsFactors=FALSE)
    if (report == "Table") Result = data.frame() else Result = vector()
    for (i in seq_len(nSUBJID)) {
      for (j in seq_len(nTRT)) {
        cSUBJID = SUBJIDs[i]
        cTRT = TRTs[j]
        cCell = (i - 1)*nTRT + j          # 1-based cell index over the SUBJID x TRT grid
        Dat = concData[concData[,id]==cSUBJID & concData[,trt]==cTRT,]
        if (nrow(Dat) > 0) {
          x = as.numeric(Dat[,Time])
          y = as.numeric(Dat[,conc])
          cDose = pick(dose, cCell)
          cTimeInfusion = pick(dur, cCell)
          if (adm == "Infusion" & !(cTimeInfusion > 0)) stop("Infusion mode should have dur larger than 0!")

          Res0 = rbind(Res0, data.frame(cSUBJID, cTRT, stringsAsFactors=FALSE))
          cRes = sNCA(x, y, dose=cDose, adm=adm, dur=cTimeInfusion, iAUC=iAUC, doseUnit=uDose, timeUnit=uTime, concUnit=uConc, down=fit, returnNA=TRUE)
          cResult = cRes[ColName0]
          if (report == "Table") {
            Result = rbind(Result, cResult)
          } else {
            Result = c(Result, "NCA REPORT", paste0("Subject=", cSUBJID), paste0("Treatment=", cTRT), cResult, "", "")
          }
        }
      }
    }
    if (report == "Table") {
      Result = cbind(Res0, Result)
      colnames(Result) = c(id, trt, ColName0)
      rownames(Result) = NULL
    }
  }

  return(Result)
}
