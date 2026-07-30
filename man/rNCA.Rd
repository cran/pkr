\name{rNCA}
\alias{rNCA}
\title{Do NCA for review}
\description{This performs NCA from the CDISC EX and PC datasets.}
\usage{
rNCA(ex, pc, study = "", trt = "", id = "", analyte = "",
     codeBQL = c("< 0", "<0", "NQ", "BLQ", "BQL", "BQoL", "<LOQ"),
     fit="Linear", MinPoints = 5, excludeInfusion = TRUE)
}
\arguments{
  \item{ex}{EX domain data, usually from the \code{loadEXPC}}
  \item{pc}{PC domain data, usually form the \code{loadEXPC}}
  \item{study}{vector of study names in EX and PC domain to do NCA. \code{""} means all.}
  \item{trt}{vector of treatment names in EXTRT to do NCA. \code{""} means all. Placebo is always excluded, since a dose of zero has no pharmacokinetics to characterize.}
  \item{id}{vector of subject IDs in USUBJID to do NCA. \code{""} means all.}
  \item{analyte}{vector of molecular species in PCTESTCD to do NCA. \code{""} means all. Matching is case-insensitive.}
  \item{codeBQL}{symbols of below the quantitation limit}
  \item{fit}{either of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC and AUMC}
  \item{MinPoints}{an analysis is done only when the number of quantifiable (greater than zero) concentrations in the dosing interval is strictly greater than this value. Note that this counts quantifiable concentrations, not sampling points: BQL and zero records do not count.}
  \item{excludeInfusion}{estimate the terminal slope of an infusion from post-infusion samples only; \code{TRUE} by default. Set \code{FALSE} for the unrestricted automatic search. See Details. No effect on a bolus or an extravascular dose.}
}
\details{
One row of result is produced for each subject, dosing interval and analyte.
A dosing interval starts at its own EXSTDTC and ends at the subject's next real
dose of any treatment in the same study, so that in a crossover trial the first
period is not extended into the second; a placebo or zero dose delivers no drug
and does not end an interval. The final interval is open-ended, so terminal and
follow-up samples still contribute to the elimination phase.

The last concentration in the 24 hours before the dose is carried into the
interval as the pre-dose sample, and placed at time zero. The lookback is bounded
because an unbounded one would report a sample from an earlier period as this
dose's initial concentration. In an interval that is followed by another dose the
trailing trough is dropped when its nominal time (PCTPTNUM) equals the pre-dose
one, so that it is not counted twice.

EX and PC records are sorted by date and time before use, so the result does not
depend on the row order of the input. Records whose date and time cannot be read
as a full \code{YYYY-MM-DDTHH:MM:SS} value are dropped with a warning naming the
subjects affected.

The automatic terminal-slope search restricts itself to samples after the maximum
concentration. That is the right rule for an extravascular dose but says nothing
about an infusion, so a sample drawn while drug is still being infused can win on
adjusted R-squared and be reported as terminal elimination. Because a terminal
phase does not exist while drug is still going in, \code{excludeInfusion} defaults
to \code{TRUE}: the search is given only the samples after the end of infusion,
and the same maximum-adjusted-R-squared criterion then selects the window among
them; within that window it still starts after the largest of those samples. Set
\code{excludeInfusion=FALSE} for the unrestricted search over all samples.

Only lambda z and the parameters derived from it (half-life, AUCIFO, AUCPEO, Vz)
are affected - Cmax, AUClast and the other observed parameters cannot change.
When there are too few post-infusion samples to fit a slope, a warning names the
profile and the unrestricted selection is used, so no profile is lost silently.

The restriction applies to this function and to \code{\link{NCA0}}, which know the
dosing record and so the end of infusion. \code{\link{sNCA}} and
\code{\link{tblNCA}} take times and concentrations without a dosing record and
remain faithful wrappers around \pkg{NonCompart}; give \code{sNCA} an explicit
\code{UsePoints} to fix the window there.

This calls NCA0.
Results of this can be further processed by foreNCA to plot and compare between studies and dose groups.
}
\value{
This returns a table of NCA results, or invisibly \code{NULL} with a warning when
no subject meets the \code{MinPoints} criterion.
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{help}}, \code{\link{NCA0}}, \code{\link{loadEXPC}}, \code{\link{foreNCA}}
}
\keyword{NCA}
