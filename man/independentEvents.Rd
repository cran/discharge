\name{independentEvents}
\alias{independentEvents}
\title{ Independent high- or low-flow events }
\description{
  Finds independent events greater than or less than a specified criterion.  High (or low) -flow days occurring on consecutive days are considered part of one event.  This function can be used to find events exceeding 2- or 10- year return levels (as calculated in \code{\link{lp3Events}} function, for example), or to find residual flows of a certain magnitude.    
}
\usage{
independentEvents(cutoff.val, data, data.column, below.cutoff=FALSE)

}
\arguments{
\item{cutoff.val}{ Numeric specifying event criterion.  }
\item{data}{ Data matrix or data frame with one column of streamflow data. }
\item{data.column}{ Numeric; specifies column in which to look for events.}
\item{below.cutoff}{ Logical. TRUE to find events less than the \code{cutoff.val} and FALSE to find events greater than the \code{cutoff.val}. }

}
\value{
A \code{data.frame} with columns
\item{events.starts}{ Index of event start. }
\item{events.ends}{ Index of event end. }
\item{events.duration}{ Length (days) of event. }
\item{extreme.this.events}{ Maximum or minimum flow for this event. }
\item{ind.extreme}{ Index of maximum or minimum flow for this event.  
If extreme is not unique, the chronologically first index is given. }
\item{...}{ All columns of original data, corresponding to max or min flow. These columns will have the same column names as the original data.}
\item{duplicates}{ 0 if the extreme is unique, 1 if it is not unique. }
}


\examples{
data(sycamore)
syc.sf<-asStreamflow(sycamore)
#find 10-year flood
q10<-lp3Events(syc.sf)$Q10

#find all events greater than 10-year flood
independentEvents(q10,syc.sf$data, data.col=8 , below.cutoff=FALSE)
}
\keyword{datasets}

\seealso{
\code{\link{lp3Events}}
}
