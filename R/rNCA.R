rNCA = function(ex, pc, study="", trt="", id="", analyte="", codeBQL=c("< 0", "<0", "NQ", "BLQ", "BQL", "BQoL", "<LOQ"), fit="Linear", MinPoints=5, excludeInfusion=TRUE)
{
# Author: Kyun-Seop Bae k@acr.kr
# Noncompartmental analysis driven by the CDISC SDTM EX and PC domains.
# One row of output per subject x dosing interval x analyte.
#
# The six near-identical code paths this function used to carry - single versus
# multiple dosing, interior versus last dosing interval, one versus several
# analytes - differed only in the row window, whether the pre-dose sample was
# recovered and whether the trailing trough was trimmed. They are now one loop
# over the dosing intervals built by .rncaIntervals(); a single dose is simply a
# subject whose only interval is also its last.
#
# INPUT
#   ex, pc: EX and PC domain tables, usually from loadEXPC()
#   study, trt, id, analyte: selections; "" means all
#   codeBQL: PCSTRESC values meaning below the limit of quantitation
#   fit: "Linear" or "Log" for the downward trapezoidal rule
#   MinPoints: an analysis needs MORE than this many quantifiable concentrations.
#     The strict comparison and the counting of quantifiable rather than sampled
#     points are both intentional and confirmed; man/rNCA.Rd documents them.
#   excludeInfusion: for an infusion, estimate the terminal slope from
#     post-infusion samples only. On by default: a terminal phase does not exist
#     while drug is still going in, so a sample drawn during the infusion is not a
#     candidate for it. Set FALSE for the unrestricted search.
# RETURNS
#   data.frame of NCA results, or invisibly NULL when nothing qualifies

  ex = .rncaDomain(ex, "ex", required=c("STUDYID", "USUBJID", "EXTRT", "EXDOSE", "EXSTDTC"),
                   defaults=list(EXDOSU="", EXROUTE="ORAL", EXENDTC=""),
                   blankAs=list(EXROUTE="ORAL"))
  pc = .rncaDomain(pc, "pc", required=c("STUDYID", "USUBJID", "PCTESTCD", "PCDTC", "PCSTRESN"),
                   defaults=list(PCSTRESU="", PCSPEC="PLASMA", PCLLOQ="", PCSTRESC=""),
                   blankAs=list(PCSPEC="PLASMA"))

  ex[toupper(ex[,"EXTRT"]) == "PLACEBO", "EXDOSE"] = "0"
# An empty EXENDTC means the dose ended when it started (mirrors readEX), so such
# records survive the datetime filter below.
  emptyEnd = !nzchar(ex[,"EXENDTC"])
  ex[emptyEnd, "EXENDTC"] = ex[emptyEnd, "EXSTDTC"]

  pc[toupper(pc[,"PCSTRESC"]) %in% toupper(Trim(codeBQL)), "PCSTRESN"] = "0"
  pc = pc[nzchar(pc[,"PCSTRESN"]), , drop=FALSE]

# Below the limit of quantitation counts as zero. Vectorised: the old row-by-row
# loop called as.numeric() twice per record and emitted a coercion warning for
# every non-numeric cell it met.
  conc = suppressWarnings(as.numeric(pc[,"PCSTRESN"]))
  lloq = suppressWarnings(as.numeric(pc[,"PCLLOQ"]))
# Reported before the rule below can turn a negative into a zero. An assay
# artefact that survives - because no PCLLOQ was given - makes every parameter of
# that profile NA, which is indistinguishable from a genuine failure to fit.
  isNeg = !is.na(conc) & conc < 0
  if (any(isNeg)) {
    warning(paste0(sum(isNeg), " PC record(s) have a negative PCSTRESN, in subject(s): ",
                   paste(utils::head(sort(unique(pc[isNeg, "USUBJID"])), 10), collapse=", ")))
  }
  pc[!is.na(conc) & !is.na(lloq) & conc < lloq, "PCSTRESN"] = "0"
  badConc = is.na(suppressWarnings(as.numeric(pc[,"PCSTRESN"])))
  if (any(badConc)) {
# One unparseable concentration used to make the MinPoints sum NA and abort the
# entire call with "missing value where TRUE/FALSE needed".
    warning(paste0(sum(badConc), " PC record(s) have a non-numeric PCSTRESN and are ignored, in subject(s): ",
                   paste(utils::head(sort(unique(pc[badConc, "USUBJID"])), 10), collapse=", ")))
  }

  study   = .rncaSelect(study, sort(unique(pc[,"STUDYID"])), "STUDYID")
  trtAll  = sort(unique(ex[,"EXTRT"]))
  trt     = .rncaSelect(trt, trtAll, "EXTRT")
  id      = .rncaSelect(id, sort(unique(pc[,"USUBJID"])), "USUBJID")
  analyte = .rncaSelect(analyte, sort(unique(pc[,"PCTESTCD"])), "PCTESTCD")
# Case variants of one code are one analyte. Deriving the list with plain
# unique() made "DRUGX" and "DrugX" two analytes and split a single profile in
# half; with the case-insensitive matching below they would instead be analysed
# twice over.
  analyte = analyte[!duplicated(toupper(analyte))]
  nAnal   = length(analyte)

# Placebo carries no drug, so there is nothing to characterise; it is excluded
# whether or not the caller named it. Saying so beats the bare "no subject" stop.
  trtU = setdiff(toupper(trt), "PLACEBO")
  if (length(trtU) == 0) stop("Only placebo was selected; placebo has no pharmacokinetics to analyze.")

  okEX = .rncaValidDTC(ex[,"EXSTDTC"]) & .rncaValidDTC(ex[,"EXENDTC"])
  .rncaWarnDropped(okEX, ex[,"USUBJID"], "EX")
  EX1 = ex[which(okEX & ex[,"STUDYID"] %in% study & toupper(ex[,"EXTRT"]) %in% trtU &
                 ex[,"USUBJID"] %in% id), , drop=FALSE]
  if (nrow(EX1) == 0) stop("No dosing record matches the given condition.")

# The route check runs on the SELECTED records and on normalised values. It used
# to run on the raw ex, so an intravenous record for a subject that would never be
# analysed - or merely the spellings "INTRAVENOUS" and "Intravenous" in one file -
# aborted the whole call.
# A blank EXROUTE was normalised to ORAL above, so it counts as extravascular
# here exactly as NCA0 will treat it. Dropping blanks before this test let an
# intravenous record and a blank-route record sit in one table as Infusion and
# Extravascular rows.
  Routes = unique(toupper(Trim(EX1[,"EXROUTE"])))
  if (length(Routes) > 1 && "INTRAVENOUS" %in% Routes) stop("Intravenous and extravascular routes cannot be mixed!")

  id = sort(unique(EX1[,"USUBJID"]))
  okPC = .rncaValidDTC(pc[,"PCDTC"])
  .rncaWarnDropped(okPC, pc[,"USUBJID"], "PC")
  PC1 = pc[which(okPC & pc[,"STUDYID"] %in% study & pc[,"USUBJID"] %in% id &
                 toupper(pc[,"PCTESTCD"]) %in% toupper(analyte) &
                 toupper(substr(pc[,"PCSPEC"], 1, 6)) == "PLASMA"), , drop=FALSE]

  IDs = sort(unique(PC1[,"USUBJID"]))
  nID = length(IDs)
  if (nID == 0) stop("No subject matches the given condition.")
  noPC = setdiff(id, IDs)
  if (length(noPC) > 0) {
    warning(paste0(length(noPC), " dosed subject(s) have no usable plasma record and are skipped: ",
                   paste(utils::head(noPC, 10), collapse=", ")))
  }

  colOut = c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")
# Doses used only to CLOSE an interval. A dose of another treatment ends the
# preceding period even when it was not selected for analysis - that is what
# keeps a crossover period from running into the next one - but it must be a real
# dose of the same study. A placebo or zero dose delivers no drug and so ends
# nothing: letting it close an interval cut AUClast by a third in a double-dummy
# design. Restricting to the selected study stops an unrelated record that merely
# shares a USUBJID from truncating, or entirely erasing, the analysis.
  exDose = suppressWarnings(as.numeric(ex[,"EXDOSE"]))
  isPlacebo = toupper(ex[,"EXTRT"]) == "PLACEBO" | (!is.na(exDose) & exDose <= 0)
  exAll = ex[which(okEX & !isPlacebo & ex[,"STUDYID"] %in% study),
             c("USUBJID", "EXSTDTC"), drop=FALSE]

# Duplicate records make the analysis ambiguous rather than wrong; say so once.
  dupEX = duplicated(paste(EX1[,"USUBJID"], EX1[,"EXSTDTC"]))
  if (any(dupEX)) {
    warning(paste0(sum(dupEX), " duplicate dosing time(s) found, in subject(s): ",
                   paste(utils::head(sort(unique(EX1[dupEX, "USUBJID"])), 10), collapse=", ")))
  }
  dupPC = duplicated(paste(PC1[,"USUBJID"], PC1[,"PCTESTCD"], PC1[,"PCDTC"]))
  if (any(dupPC)) {
    warning(paste0(sum(dupPC), " duplicate sampling time(s) found, in subject(s): ",
                   paste(utils::head(sort(unique(PC1[dupPC, "USUBJID"])), 10), collapse=", ")))
  }
  Rows = list()
  for (i in seq_len(nID)) {
    cID = IDs[i]
    EXi = unique(EX1[EX1[,"USUBJID"] == cID,
                     c("STUDYID", "USUBJID", "EXTRT", "EXDOSE", "EXDOSU", "EXROUTE", "EXSTDTC", "EXENDTC")])
    if (nrow(EXi) == 0) next
    Ivls = .rncaIntervals(EXi, exAll[exAll[,"USUBJID"] == cID, "EXSTDTC"])
    PCi = PC1[PC1[,"USUBJID"] == cID, , drop=FALSE]

    for (Ivl in Ivls) {
      cEX = Ivl$ex
      for (k in seq_len(nAnal)) {
# Case-insensitive, matching the PCTESTCD filter that built PC1. The inner loops
# used to compare exactly, so analyte="druga" against data holding "DRUGA" passed
# the outer filter and then selected nothing at all.
        cAnalyte = analyte[k]
        PCk = PCi[which(toupper(PCi[,"PCTESTCD"]) == toupper(cAnalyte)), , drop=FALSE]
        if (nrow(PCk) == 0) next
# Sort by time. NonCompart::sNCA rejects unsorted times outright, and the
# pre-dose and trough rules below both read the first and last record. The
# concentration breaks ties, so two records sharing a sampling time cannot make
# the result depend on which of them came first in the file.
        PCk = PCk[order(PCk[,"PCDTC"], suppressWarnings(as.numeric(PCk[,"PCSTRESN"]))), , drop=FALSE]

        PCy = .rncaWindow(PCk, Ivl)
        PCy = .rncaTrimTrough(PCy, Ivl)
        if (nrow(PCy) == 0) next
        if (sum(suppressWarnings(as.numeric(PCy[,"PCSTRESN"])) > 0, na.rm=TRUE) <= MinPoints) next

# NCA0 indexes its EX argument by name, so it needs a named character vector.
# Every call site used to hand it a one-row data.frame, whose [ ] returns a
# one-column data.frame rather than a string; that reached NonCompart::Unit() and
# failed with "non-character argument". It is the reason rNCA could not run at
# all on the data.frames loadEXPC() produces.
        Res0 = NCA0(unlist(lapply(cEX[1,], as.character)), PCy[, colOut, drop=FALSE], fit=fit,
                    excludeInfusion=excludeInfusion,
                    label=paste0(cID, " / ", cAnalyte, " / dose at ", Ivl$start))
# The label comes from the data, not from the caller's spelling, so that
# downstream matching on PCTESTCD (foreNCA) finds the rows it expects.
        cLabel = unique(PCy[,"PCTESTCD"])[1]
        Res1 = c(cEX[1,"STUDYID"], cEX[1,"USUBJID"], cEX[1,"EXTRT"], cLabel, Ivl$start, Res0)
        names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", names(Res0))
        Rows[[length(Rows) + 1]] = Res1
      }
    }
  }

  if (length(Rows) == 0) {
    warning("No subject passed the minimum-points criterion; nothing to analyze.")
    return(invisible(NULL))
  }

  Res = .rncaBindRows(Rows)
  Res = Res[order(Res[,"STUDYID"], Res[,"USUBJID"], Res[,"PCRFTDTC"], Res[,"PCTESTCD"]), , drop=FALSE]
  Res = conv.pp(Res)
  attr(Res, "NCA") = "ncaRes"
  class(Res) = union(class(Res), "pp")
  return(Res)
}
