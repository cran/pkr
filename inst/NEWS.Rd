\name{NEWS}
\title{News for Package \pkg{pkr}}
\section{Version 0.1.6 (2026-07-30)}{
  \itemize{
    \item{\code{plotPK} no longer loses series when a study has more than 18 subjects or treatments. The colour vector had a fixed length of 18, and indexing past its end gave \code{NA}, which draws nothing at all, so every subject after the 18th was silently missing from the pooled profiles. Colours are now generated for any number of series, and plotting symbols vary as well.}
    \item{\code{plotPK} draws the subject legend in a panel of its own, sizing the text and the number of columns to fit. Previously the legend was placed inside the plot region with a fixed colour vector, so with many subjects it overlapped the curves, ran off the figure, and mislabelled series through silent recycling.}
    \item{\code{plotPK} lays treatment panels out in a grid of at most four columns. The former single row of \code{nTRT} panels asked for a right margin of \code{2.5*nTRT} lines and a raster \code{700*nTRT} pixels wide, which failed with \dQuote{figure margins too large} from about eight treatments onward.}
    \item{\code{plotPK} drops non-finite values before computing axis ranges. A single \code{NA} concentration used to abort the whole function with \dQuote{need finite 'ylim' values}.}
    \item{\code{plotPK} closes every graphics device it opens, even when a panel fails. A failed call used to leave the device open, so the caller's next plot was written into the half-finished file.}
    \item{\code{plotPK} sanitizes the file-name stem taken from the \code{concData} argument, and accepts an explicit \code{name} argument. Calling it with an expression such as \code{d[d$ID != "S001", ]} used to fail with \dQuote{cannot open file}.}
    \item{\code{plotPK} prints the per-bin counts as one line; the whole vector used to be drawn on top of itself. Its mean profile also sizes the y axis to the error bars instead of clipping them, warns about subjects left out of the mean, and no longer requires \pkg{binr}.}
    \item{Semi-logarithmic axes in \code{plotPK}, \code{pdfNCA} and \code{rtfNCA} are aligned to whole decades enclosing the data. The previous \code{trunc()}-based tick positions could all fall outside \code{ylim}, leaving the y axis blank.}
    \item{\code{plotFit} plots only positive concentrations instead of substituting 0.1 for zeros, includes the fitted intercept in the y range, and reports a clear error for an unknown \code{ID}.}
    \item{\code{conv.pp} works with a single-column result; it used to fail with \dQuote{subscript out of bounds}.}
    \item{\code{Subset} compares whole columns instead of one cell at a time, and rejects an unknown key with a clear message.}
    \item{\code{foreNCA} checks for the columns it needs and reports an empty selection instead of failing inside \pkg{forestplot}.}
    \item{\code{NCA0} uses a single concentration unit and warns when the PC records disagree. \code{readEX} no longer requires an \code{EXENDTC} column.}
    \item{\code{plotPK} returns the paths it wrote, invisibly.}
    \item{\code{rNCA} now runs at all. Since the first release it has passed a one-row data.frame to \code{NCA0}, which indexes its EX argument by name; the result was a one-column data.frame where a string was expected, and every call made the way the documentation prescribes - with \code{ex} and \code{pc} from \code{loadEXPC()} - died inside \code{NonCompart::Unit()} with \dQuote{non-character argument}. Passing a matrix instead failed differently, because subsetting a matrix drops a single matching row to a vector. Both input shapes are now accepted and normalized once on entry.}
    \item{\code{rNCA}'s six near-duplicate code paths (single versus multiple dosing, interior versus last dosing interval, one versus several analytes) are one loop over dosing intervals. A single dose is simply a subject whose only interval is also its last.}
    \item{\code{rNCA} sorts EX and PC records by date and time, so the result no longer depends on the row order of the input. Previously \code{unique()} preserved file order, so EX rows supplied in another order - which \code{combXPT} can produce, since it binds folders in directory order - silently merged two dosing intervals into one and dropped the other, and unsorted PC rows either produced a scrambled profile or failed inside \pkg{NonCompart}.}
    \item{\code{rNCA} ends each dosing interval at the subject's next dose of any treatment. In a crossover trial analysed one treatment at a time, the first period used to run open-ended and absorb the second period's samples, overstating AUClast about threefold and taking Cmax from the wrong period.}
    \item{\code{rNCA} binds result rows by name. \code{rbind} matches by position and takes its column names from the first row, so a cohort mixing administration methods - whose parameter sets differ in length - wrote later rows' values under the wrong headers, reporting for instance a volume of distribution as a clearance, with only a recycling warning.}
    \item{\code{rNCA} finds the pre-dose sample as the last record in the day before the dose rather than by matching calendar dates. The old rule made the result depend on the wall-clock time of dosing: a dose shortly after midnight could never recover a trough drawn the previous evening. The lookback is bounded in real time because \code{NCA0} places a pre-dose record at time zero, so an unbounded search would report a sample from the previous period as this dose's initial concentration.}
    \item{\code{rNCA} closes a dosing interval only at a real dose of the study being analysed. A placebo or zero dose delivers no drug and no longer truncates an interval, and a record that merely shares a subject identifier with another study no longer does either.}
    \item{\code{rNCA} treats a blank EXROUTE as extravascular, as \code{NCA0} always has, so a cohort mixing intravenous with blank-route records is refused rather than silently reported with two different models side by side. A blank PCSPEC likewise now means plasma, as an absent PCSPEC column already did.}
    \item{\code{rNCA} orders duplicate sampling and dosing times deterministically and reports them, so a repeated PCDTC or EXSTDTC can no longer make the result depend on the row order of the input. Negative concentrations are also reported.}
    \item{\code{rNCA} resolves \code{study}, \code{trt}, \code{id} and \code{analyte} to the spellings found in the data, so every selector matches case-insensitively end to end. \code{id=c("SUBJ-01","subj-02")} used to return only the first subject, silently.}
    \item{\code{rNCA} no longer drops the last sample of every interior dosing interval when a PCTPTNUM column is present but blank, which the old unguarded string comparison did silently.}
    \item{\code{rNCA} matches analyte names case-insensitively throughout, and labels the output PCTESTCD from the data rather than from the caller's spelling. \code{analyte="druga"} against data holding \code{"DRUGA"} used to select nothing, and case variants of one code were analysed as two half profiles.}
    \item{\code{rNCA} tolerates a non-numeric PCSTRESN, warning and ignoring it, instead of aborting the whole run with \dQuote{missing value where TRUE/FALSE needed}. It also accepts \code{character(0)} and \code{NULL} selections, checks for the SDTM columns it needs, supplies the same defaults as \code{readEX}/\code{readPC} for the optional ones, warns when records are dropped for an unusable date/time or when a dosed subject has no plasma record, and checks routes after selection rather than before.}
    \item{\code{rNCA} and \code{NCA0} estimate the terminal slope of an infusion from post-infusion samples only, controlled by the new \code{excludeInfusion} argument and ON BY DEFAULT. The automatic search constrains itself to samples after the maximum concentration, which says nothing about an infusion, so a sample drawn while drug was still going in could win on adjusted R-squared and be reported as terminal elimination; a terminal phase does not exist while drug is still going in. The same maximum-adjusted-R-squared criterion still selects the window, now among the post-infusion samples. Only lambda z and what derives from it (half-life, AUCIFO, AUCPEO, Vz) are affected; Cmax and AUClast cannot change. Too few post-infusion samples produces a warning naming the profile and falls back to the unrestricted search, so no profile is lost silently. Pass \code{excludeInfusion=FALSE} to reproduce results from an earlier version. The restriction applies to the SDTM entry points, which know the dosing record; \code{sNCA} and \code{tblNCA} remain faithful wrappers around \pkg{NonCompart}.}
    \item{\code{sNCA} gained \code{UsePoints}, passed through to \pkg{NonCompart}, which fixes the terminal-slope window explicitly instead of searching for it. The \code{UsedPoints} attribute now survives the \code{returnNA=TRUE} padding.}
    \item{\code{NCA0} decides the administration method from EXROUTE first. Testing the recorded duration first classified any extravascular record carrying an EXENDTC - routine in real EX domains - as an infusion, so an oral dose was reported with VZO/CLO/VSS in place of VZFO/CLFO/MRTEV.}
  }
}
\section{Version 0.1.5 (2026-07-22)}{
  \itemize{
    \item{Core noncompartmental computation is now delegated to the \pkg{NonCompart} package, so \pkg{pkr} and \pkg{NonCompart} share one maintained engine and always give identical results. This fixes, among others, the handling of \code{NA} concentrations in \code{sNCA} and the column mis-alignment of \code{tblNCA}/\code{NCA} with heterogeneous subjects.}
    \item{\code{NCA()} with a treatment (\code{trt}) argument no longer errors on an undefined dose object.}
    \item{\code{readPC} now defaults optional \code{PCSPEC}/\code{PCLLOQ}/\code{PCSTAT}/\code{PCORRES} columns correctly and recognizes more BLQ codes.}
    \item{\code{combXPT} no longer requires the first folder to contain the domain file.}
    \item{\code{plotPK} now writes to the \code{outdir} argument, and its mean annotations use the actual means.}
    \item{\code{txtNCA} reports the correct administration route and infusion length.}
    \item{\code{plotFit} draws the semi-log axis and terminal-slope line correctly.}
    \item{\pkg{foreign}, \pkg{binr}, \pkg{forestplot}, and \pkg{rtf} moved from Depends to Imports; \pkg{NonCompart} added to Imports.}
  }
}
\section{Version 0.1.4 (2022-07-10)}{
  \itemize{
    \item{Rare cases raising NaN in \code{BestSlope} are handled.}
  }
}
\section{Version 0.1.3 (2020-06-10)}{
  \itemize{
    \item{Column order bug with bolus infusion fixed.}
  }
}
\section{Version 0.1.2 (2018-06-04)}{
  \itemize{
    \item{Deprecated option in \code{Sys.timezone()} is removed.}
  }
}
\section{Version 0.1.1 (2017-07-10)}{
  \itemize{
    \item{\pkg{pkr} A bug in C0 calcuation of IV bolus case was fixed.}
  }
}
\section{Version 0.1.0 (2017-03-15)}{
  \itemize{
    \item{First version}
  }
}
