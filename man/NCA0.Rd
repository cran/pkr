\name{NCA0}
\alias{NCA0}
\title{NCA of SDTM data for single subject}
\description{
This performs Noncompartmental Analysis(NCA) for one dosing record from the CDISC EX and PC domain.
}
\usage{NCA0(EX0, PC0, fit="Linear", excludeInfusion=TRUE, label="")}
\arguments{
  \item{EX0}{One EX record: a named character vector, a one-row data.frame or a one-row matrix. It needs the names EXDOSE, EXDOSU, EXROUTE, EXSTDTC and EXENDTC.}
  \item{PC0}{The concentration records of that dosing interval, with the columns PCDTC, PCSTRESN and PCSTRESU}
  \item{fit}{either of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC and AUMC}
  \item{excludeInfusion}{when the administration is an infusion, estimate the terminal slope from post-infusion samples only; \code{TRUE} by default. See \code{\link{rNCA}} for what this changes.}
  \item{label}{identifier used only in warning messages}
}
\details{
The administration method is taken from EXROUTE first: an EXROUTE of
\code{"INTRAVENOUS"} gives \code{"Infusion"} when EXENDTC is later than EXSTDTC
and \code{"Bolus"} otherwise, and any other route gives \code{"Extravascular"}.
Time after dose is computed from EXSTDTC and negative values are set to zero, so
a pre-dose sample is placed at time zero.

This calls \code{sNCA}.
This is called by \code{rNCA} function.
}
\value{
This returns a named character vector: AdmMethod, Dose, UnitDose, Dur, UnitTime,
UnitConc, followed by the NCA parameters. Which parameters are present depends on
the administration method.
}
\author{Kyun-Seop Bae <k@acr.kr>}

\seealso{
\code{\link{help}}, \code{\link{rNCA}}, \code{\link{sNCA}}
}
\keyword{ NCA }
