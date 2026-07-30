# Author: Kyun-Seop Bae k@acr.kr
# Internal helpers for the SDTM-driven analysis in rNCA(). Their names start with
# a dot, which the exportPattern in NAMESPACE does not match, so they stay out of
# the public API.

.rncaDomain = function(d, what, required, defaults=list(), blankAs=list())
{
# Normalise one SDTM domain: a data.frame of trimmed character columns with NA
# rendered as "". rNCA used to accept whatever it was handed, but a matrix and a
# data.frame behave differently under [ , ] (a matrix drops to a vector when one
# row matches, a data.frame does not), and NCA0 needs a named character vector.
# Coercing once here is what lets the rest of the code have a single input shape.
  if (is.null(d) || length(dim(d)) != 2) stop(paste0(what, " should be a data.frame or a matrix!"))
  d = as.data.frame(d, stringsAsFactors=FALSE)
  if (nrow(d) == 0) stop(paste0(what, " has no record!"))

  missingCol = setdiff(required, colnames(d))
  if (length(missingCol) > 0) {
    stop(paste0(what, " has no column: ", paste(missingCol, collapse=", "), "!"))
  }

  for (i in seq_len(ncol(d))) {
    d[,i] = Trim(as.character(d[,i]))
    d[is.na(d[,i]), i] = ""
  }
# Optional columns get the same defaults readEX()/readPC() would have supplied, so
# rNCA behaves identically whether or not the caller came through loadEXPC().
  for (cName in names(defaults)) {
    if (!(cName %in% colnames(d))) d[,cName] = defaults[[cName]]
  }
# A column that is present but blank must mean the same as an absent one.
# Otherwise a blank EXROUTE slipped past the route-mixing check while NCA0 still
# modelled it as extravascular, and a blank PCSPEC failed the PLASMA test and
# emptied the analysis with a message blaming the caller's selection.
  for (cName in names(blankAs)) {
    if (cName %in% colnames(d)) d[!nzchar(d[,cName]), cName] = blankAs[[cName]]
  }
  return(d)
}

.rncaSelect = function(arg, all, what)
{
# Resolve a selection argument. "" (the documented "everything" marker), NULL and
# character(0) all mean "all"; the old `if (arg[1] == "")` raised "missing value
# where TRUE/FALSE needed" for the latter two, which programmatic callers hit
# whenever their filter matched nothing.
  if (length(arg) == 0 || !nzchar(arg[1])) return(all)
# Return the values as they are SPELT IN THE DATA. The match is case-insensitive,
# but the filters downstream are not, so returning the caller's spelling made
# id=c("SUBJ-01","subj-02") quietly analyse only the first of the two.
  hit = toupper(all) %in% toupper(arg)
  if (!any(hit)) stop(paste0("No ", what, " in the data matches: ", paste(arg, collapse=", "), "!"))
  miss = arg[!(toupper(arg) %in% toupper(all))]
  if (length(miss) > 0) {
    warning(paste0("No ", what, " in the data matches: ", paste(miss, collapse=", "), "; ignored."))
  }
  return(all[hit])
}

.rncaValidDTC = function(x)
{
# rNCA has always kept only full 19-character ISO datetimes, because the window
# comparisons below are lexicographic and only order-correct at a fixed width.
# That filter is preserved; what is new is that dropping a record is no longer
# silent - a partial datetime used to turn a multi-dose subject into a
# single-dose one with no message at all.
  return(!is.na(x) & nchar(x) == 19 & !is.na(strptime(x, "%Y-%m-%dT%H:%M:%S")))
}

.rncaWarnDropped = function(keep, ids, what)
{
  n = sum(!keep)
  if (n > 0) {
    who = sort(unique(ids[!keep]))
    warning(paste0(n, " ", what, " record(s) dropped for an unusable date/time, in subject(s): ",
                   paste(utils::head(who, 10), collapse=", "),
                   if (length(who) > 10) ", ..." else ""))
  }
  return(invisible(NULL))
}

.rncaIntervals = function(EXi, allStart=character(0))
{
# The dosing intervals of one subject, in chronological order.
#
# EXi arrives from unique(), which preserves file order rather than time order,
# so without this sort EXi[j+1,] was not necessarily the next dose: a legal SDTM
# extract with its EX rows in another order silently produced one merged interval
# and dropped the other, with no warning. combXPT() rbinds folders in directory
# order, so that is not a hypothetical.
#
# An interval ends at the subject's next dose of ANY treatment (allStart), not
# merely the next SELECTED one. In a crossover trial analysed one treatment at a
# time, bounding by the selected records alone left the first period open-ended,
# so the second period's samples were folded into it - measured at 2.6 times the
# true AUClast, with Cmax taken from the wrong period.
# The tie-break keeps the order deterministic when a subject has two dose records
# at the same instant, which otherwise left the output row order dependent on the
# order of the input file.
  EXi = EXi[order(EXi[,"EXSTDTC"], EXi[,"EXTRT"], EXi[,"EXDOSE"]), , drop=FALSE]
  allStart = sort(unique(c(allStart, EXi[,"EXSTDTC"])))
  Out = vector("list", nrow(EXi))
  for (j in seq_len(nrow(EXi))) {
    tStart = EXi[j, "EXSTDTC"]
    Next = allStart[allStart > tStart]
    Prev = allStart[allStart < tStart]
    Out[[j]] = list(ex    = EXi[j, , drop=FALSE],
                    start = tStart,
# The final interval stays open-ended, as it always has been, so that terminal
# and follow-up samples still contribute to the elimination phase.
                    end   = if (length(Next) > 0) Next[1] else NA_character_,
                    prev  = if (length(Prev) > 0) Prev[length(Prev)] else NA_character_,
                    last  = (length(Next) == 0))
  }
  return(Out)
}

.rncaWindow = function(PCi, Ivl, predoseHours=24)
{
# The concentration records belonging to one dosing interval: everything strictly
# after this dose and, for an interior interval, strictly before the next one -
# plus the pre-dose sample.
#
# The pre-dose sample is the last record in the predoseHours before the dose.
# rNCA used to look for it among the records sharing the dose's CALENDAR DATE,
# which made the result depend on the wall-clock time of dosing: a dose at 00:15
# could never recover a trough drawn the previous evening, so two identical
# designs differing only in clock time reported different AUCs.
#
# The lookback has to be bounded in real time, not merely by the previous dose.
# NCA0 clamps a negative time after dose to zero, so any record admitted here is
# injected as a concentration at time zero. With only the previous dose as the
# bound, a crossover period whose own pre-dose sample was not drawn picked up the
# LAST sample of the preceding period - 70 h earlier, near that period's peak -
# and reported it as this period's C(0): Cmax 7.12 at Tmax 0 against a true 4.75
# at 6 h, and AUClast 48 percent high. One day matches the intent of the original
# calendar-date rule without inheriting its dependence on the clock, and is the
# bound the maintainer settled on; do not widen it without revisiting that case.
  after = PCi[,"PCDTC"] > Ivl$start
  if (!is.na(Ivl$end)) after = after & PCi[,"PCDTC"] < Ivl$end
  PCy = PCi[which(after), , drop=FALSE]

  Floor = format(as.POSIXct(Ivl$start, format="%Y-%m-%dT%H:%M:%S", tz="UTC") - predoseHours*3600,
                 "%Y-%m-%dT%H:%M:%S", tz="UTC")
  cand = PCi[,"PCDTC"] <= Ivl$start & PCi[,"PCDTC"] >= Floor
  if (!is.na(Ivl$prev)) cand = cand & PCi[,"PCDTC"] > Ivl$prev
  Pre = PCi[which(cand), , drop=FALSE]
  if (nrow(Pre) > 0) PCy = rbind(Pre[nrow(Pre), , drop=FALSE], PCy)

  return(PCy)
}

.rncaTrimTrough = function(PCy, Ivl)
{
# In an interior interval the final record is the next dose's pre-dose trough,
# which carries the same nominal time as this interval's own pre-dose sample;
# rNCA has always dropped it so that the trough is not counted twice. Only the
# last interval, which has no following dose, keeps its final record.
#
# The comparison is now guarded. It used to be a bare string equality, so a
# PCTPTNUM column that was present but blank made "" == "" true and silently
# removed the last sample of EVERY interior interval - measured at 24% low on
# AUClast. Blank labels no longer match, and "0" and "0.0" now do.
  if (isTRUE(Ivl$last)) return(PCy)
  if (!("PCTPTNUM" %in% colnames(PCy)) || nrow(PCy) < 2) return(PCy)
  a = PCy[1, "PCTPTNUM"]
  b = PCy[nrow(PCy), "PCTPTNUM"]
  if (!nzchar(a) || !nzchar(b)) return(PCy)
  na = suppressWarnings(as.numeric(a))
  nb = suppressWarnings(as.numeric(b))
  same = if (is.na(na) || is.na(nb)) identical(a, b) else isTRUE(all.equal(na, nb))
  if (same) PCy = PCy[seq_len(nrow(PCy) - 1), , drop=FALSE]
  return(PCy)
}

.postInfusionPoints = function(x, y, dur, adm)
{
# Indices of the points sNCA should use for the terminal slope when the terminal
# phase is to be taken from after the end of infusion only.
#
# NonCompart::BestSlope constrains its search to points after the maximum
# concentration, which is the right rule for an extravascular profile but says
# nothing about an infusion: a sample drawn while drug is still going in can win
# on adjusted R-squared and be reported as terminal elimination. Here the search
# is handed only the post-infusion points, and BestSlope's own criterion picks the
# window inside them, so no selection rule is re-implemented.
#
# Returns indices into the NA-removed x/y that sNCA sees, suitable for its
# UsePoints argument, or NULL when no admissible window exists.
  if (!isTRUE(dur > 0)) return(NULL)
  keep = !is.na(x) & !is.na(y)
  x = x[keep]
  y = y[keep]
  cand = which(x > dur & y > 0 & is.finite(y))
# BestSlope needs at least three points to return a fitted window, and it starts
# after the maximum of whatever it is given, so a fourth is needed in practice.
  if (length(cand) < 4) return(NULL)
  bs = try(NonCompart::BestSlope(x[cand], y[cand], adm=adm), silent=TRUE)
  if (inherits(bs, "try-error") || !isTRUE(bs["LAMZNPT"] >= 2)) return(NULL)
  lo = which(x[cand] == bs["LAMZLL"])[1]
  hi = which(x[cand] == bs["LAMZUL"])[1]
  if (is.na(lo) || is.na(hi) || hi <= lo) return(NULL)
  return(cand[lo:hi])
}

.rncaBindRows = function(Rows)
{
# Bind the per-analysis result vectors BY NAME. rbind() matches by position and
# takes its column names from the first row bound, so a cohort mixing
# administration methods - whose parameter sets differ in length - wrote each
# subsequent row's values under the wrong headers. That produced, for example, a
# volume of distribution printed as a clearance, with nothing worse than a
# "number of columns of result is not a multiple of vector length" warning.
  allNames = unique(unlist(lapply(Rows, names)))
  Out = matrix(NA_character_, nrow=length(Rows), ncol=length(allNames),
               dimnames=list(NULL, allNames))
  for (i in seq_along(Rows)) {
    r = Rows[[i]]
    Out[i, names(r)] = as.character(r)
  }
  return(Out)
}
