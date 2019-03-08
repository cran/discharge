\name{allstats}
\alias{allstats}
\title{ Calculate all characteristic stats for a site }
\description{
  Calculates all parameters at once for a single data file.  The output gives numerical results from the functions \code{\link{annualnoise}}, \code{\link{fourierAnalysis}}, \code{\link{lp3Events}},\code{\link{sigmaHighFlows}}, and \code{\link{sigmaLowFlows}}.  
}
\usage{
allstats(file.name, river.name, file.type="txt", date.col=3, 
discharge.col=4, skipped.rows=28)
}
\arguments{
\item{file.name}{ Character string of the form "file.txt" or "file.csv".  }
\item{river.name}{ Character string specifying river name.}
\item{file.type}{Character string, "txt" or "csv".  Defaults to "txt".}
\item{date.col}{ Numeric specifying column containing date in "MM-DD-YYYY" format.  Defaults to 3.}
\item{discharge.col}{ Numeric specifying column containing discharge data. Defaults to 4. }
\item{skipped.rows}{Numeric indicating number of rows to skip at beginning of file.}
}
\value{
A data frame with columns
\item{a.rms}{ Root mean squared amplitude. }
\item{n.rms}{ Root mean squared noise. }
\item{snr}{ Signal-to-noise ratio. }
\item{theta.d}{ Daily noise color. }
\item{theta.a}{ Annual noise color. }
\item{sigma.lf}{ Sigma for low flow events. }
\item{sigma.hf}{ Sigma for high flow events. }
\item{q2}{ 2-year return level (flood).}
\item{q10}{ 10-year return level (flood).}
\item{l2}{ 2-year return level (drought).}
\item{l10}{ 10-year return level (drought).}
}



\note{
The arguments "date.col", "discharge.col", and "skipped.rows" are designed to give some flexibility in file input; however, tab-delimited text without extra columns will likely work best. 
}

\examples{
# allstats function works on files
# read R data into temporary file handle
data(sycamore)
f = tempfile(fileext="txt")
write.table(sycamore, file=f, sep="\t")

# print all statistics for this river
allstats(f,river.name="sycamore", date.col = 2, discharge.col = 3, skipped.rows = 1)

}

\seealso{
\code{\link{parameters.list}}
}

