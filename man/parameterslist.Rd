\name{parameters.list}
\alias{parameters.list}
\title{ Calculate all characteristic stats from a list of files }
\description{
  Function takes a vector of file names and returns seasonal signal-to-noise ratio, daily and annual noise color, 2- and 10-year return levels, and \eqn{\sigma} for high- and low-flow events.  
  
  All files in the vector should be of the same format.
}
\usage{
parameters.list(x, names=NULL, file.type="txt", date.col=3, dis.col=4, skipped.rows=28) 

}
\arguments{
\item{x}{ Character vector containing file names.}
\item{names}{ Optional character vector with names of sites.}
\item{file.type}{Character string, "txt" or "csv".  Defaults to "txt".}
\item{date.col}{ Numeric specifying column containing date in "MM-DD-YYYY" format.  Defaults to 3.}
\item{dis.col}{ Numeric specifying column containing discharge data. Defaults to 4. }
\item{skipped.rows}{Numeric indicating number of rows to skip at beginning of file.}
}
\value{
A data frame with one row for each file and the following columns:
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
The arguments "date.col", "discharge.col", and "skipped.rows" are designed to give some flexibility in file input; however, tab-delimited text without extra columns will work best. 
}
\keyword{datasets}

\examples{
# this function works on list of files
# read R data into temporary file handle
data(sycamore)
f = tempfile(fileext="txt")
write.table(sycamore, file=f, sep="\t")

# print all statistics for the list of rivers
parameters.list(c(f), names=c("sycamore"), date.col=2,dis.col=3,skipped.rows = 1)

}


\seealso{
\code{\link{allstats}}
}
