\name{plotPK}
\alias{plotPK}
\title{
Plot concentration vs. time curve for individuals and collectively.
}
\description{
Generates individual and superposed concentration vs. time curves and saves them in pdf and tiff files.
}
\usage{
plotPK(concData, id, Time, conc, unitTime = "hr", unitConc = "ng/mL", trt = "",
       fit = "Linear", dose = 0, adm = "Extravascular", dur = 0, outdir = "Output",
       name = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{concData}{name of data table containing time-concentration data of multiple subjects}
  \item{id}{column name for subject ID}
  \item{Time}{column name for the time}
  \item{conc}{column name for the concentration}
  \item{unitTime}{unit for the time}
  \item{unitConc}{unit for the concentration}
  \item{trt}{column name for the treatment code. This is useful for crossover study like bioequivalence trial.}
  \item{fit}{one of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC}
  \item{dose}{administered dose. One should be careful for the unit. This can be a vector containing dose for each subject and treatment in order.}
  \item{adm}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
  \item{dur}{infusion duration for constant infusion, otherwise 0.  This can be a vector containing values for each subject and treatment in order.}
  \item{outdir}{name of the folder to be used for the output files. \code{""} means the current working directory.}
  \item{name}{stem used in the output file names. \code{""} means the deparsed \code{concData} argument, sanitized so that it is a valid file name.}
}
\details{
This function generates five figures: individual linear and semi-logarithmic profiles
(one pdf each, several panels per page) and pooled linear, semi-logarithmic and
mean with 95\% confidence interval profiles (one tiff each, or png where the R build
has no tiff device).

Each subject gets its own colour and plotting symbol whatever the number of subjects,
and the subject legend is drawn in a panel of its own so that it cannot overlap the
curves. When there are more subjects than can be labelled legibly, the count is shown
instead of an unreadable legend. Treatment panels are laid out in a grid of at most
four columns, so a study with many treatments neither exhausts the figure margins nor
produces an enormous raster.

Non-finite times and concentrations are dropped from the axis ranges, and the
semi-logarithmic figures use only positive concentrations. The semi-logarithmic axis
is aligned to whole decades enclosing the data, so a labelled tick is always visible.

This function calls \code{NCA()}.
}
\value{
Invisibly, a character vector of the paths written. The files are saved in the
\code{outdir} folder.
}
\author{Jee Eun Lee <JeeEun.Lee@fda.hhs.gov>, Kyun-Seop Bae <k@acr.kr>}

\seealso{\code{\link{NCA}}}
\examples{
# plotPK(Theoph, "Subject", "Time", "conc", unitTime="hr", unitConc="mg/L", dose=320)
# plotPK(Indometh, "Subject", "time", "conc", unitTime="hr", unitConc="mg/L", adm="Bolus", dose=25)
}
\keyword{Plot}
