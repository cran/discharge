\name{asStreamflow}
\alias{asStreamflow}
\alias{print.streamflow}
\alias{summary.streamflow}
\title{ Create streamflow object }
\description{
  This function converts a dataset to S3 object of class "streamflow".  Data in this format can be provided as arguments to other functions to call default procedures.
}
\usage{
asStreamflow(x,river.name=NULL, start.date=NULL, end.date=NULL,max.na=10)
\method{print}{streamflow}(x, ...)
\method{summary}{streamflow}(object, ...)
}
\arguments{
\item{x,object}{ A matrix with dates in the first column and discharge values in the second column.  Dates should be of the format "YYYY-MM-DD".  }
\item{river.name}{ A character vector listing the river name. }
\item{start.date}{ Optional character string giving date to start analysis, of the format "YYYY-MM-DD"}
\item{end.date}{ Optional date to start analysis, of the format "YYYY-MM-DD" }
\item{max.na}{ Optional number specifying maximum NA values to allow. }
\item{...}{Other parameters. }
}
\value{
An object of class \code{streamflow} containing the following items:
\item{data}{ Data frame with 8 columns }
\item{n}{ Number of observations in series. }
\item{n.nas}{ Number of NA observations in series. }
\item{start}{ Start date. }
\item{end}{ End date. }
\item{name}{ Name of river. }
\item{na.info}{ Matrix containing the index and start date of all blocks of more than one NA observation. }


}

\examples{
data(sycamore)
sycamore.flows<-asStreamflow(sycamore,start.date="1965-01-01", 
                             end.date="2010-03-15",river.name="Sycamore Creek")
}
