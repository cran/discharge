\name{lp3Events}
\alias{lp3Events}
\title{ Find 2- and 10-year return levels }
\description{
  Uses the method of moments to find the 2- and 10-year droughts and floods under the Log-III Pearson distribution. 
}
\usage{
lp3Events(x)

}
\arguments{
\item{x}{ Object of class "\code{streamflow}". }
}
\value{
A \code{list} with items
\item{Q2}{ 2-year high return level. }
\item{Q10}{ 10-year high return level. }
\item{L2}{ 2-year low return level. }
\item{L10}{ 10-year low return level. }

}
\details{
Return levels are calculated using the method of moments through the package \href{https://cran.r-project.org/package=lmom}{"lmom"}.  High return levels (floods) are calculated using annual maxima of the raw (log normalized) data, while the low return levels (droughts) are calculated using annual mimima of the 7-day moving averages. 
}

\examples{
data(sycamore)
syc.sf<-asStreamflow(sycamore)
lp3Events(syc.sf)

}

