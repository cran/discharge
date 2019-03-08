\name{annualExtremes}
\alias{annualExtremes}
\title{ Calculate annual extreme flows. }
\description{
  Calculates annual extreme events for all years in the series.  By default, the function finds annual extreme discharge in \code{streamflow} object, but any matrix or \code{data.frame} may be used.
}
\usage{
annualExtremes(x,data.col=NULL, year.col=NULL, moving.avg=FALSE)

}
\arguments{
\item{x}{ Object from which to extract extremes. Should be of class \code{streamflow} or \code{data.frame} or matrix.  }
\item{data.col}{ Optional.  If input is a matrix or \code{data.frame}, specifies which column contains the data.}
\item{year.col}{ Optional.  If input is a matrix or \code{data.frame}, specifies which column contains the year.}
\item{moving.avg}{ Logical; defaults to FALSE. Can be specified TRUE to use 7-day moving average discharge when input is of class "\code{streamflow}". }
}
\value{
A \code{list} with items
\item{annual.max}{ Matrix giving maximum flow for each year in series.  Each row contains the maximum values and all corresponding variables from that observation.}
\item{annual.min}{ Matrix giving minimum flow for each year in series.  Each row contains the minimum values and all corresponding variables from that observation.}

}
\examples{
data(sycamore)
sycamore.flows<-asStreamflow(sycamore,river.name="Sycamore Creek")
syc.extremes<-annualExtremes(sycamore.flows)
names(syc.extremes)
syc.extremes$annual.max[1:3,]

}

