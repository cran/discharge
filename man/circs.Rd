\name{circ.s}
\alias{circ.s}
\title{ Estimate directional statistics for one-sigma events }
\description{
Uses the package \href{https://cran.r-project.org/package=CircStats}{"CircStats"} to find the mean direction, rho, and kappa of the von Mises distribution summarizing the ordinal day of high- and low-flow events.
}
\usage{
circ.s(x) 

}
\arguments{
\item{x}{ Output from "compare.periods" function (of class "\code{compflows}"). }
}
\value{
\item{circstats}{ \code{data.frame} with rows corresponding to the high- and low-flows for both periods. The columns list n, mu, rho, and kappa as calculated using the CircStats package. }

}

\examples{
# load data
data("sycamore")

# compare for periods from 1960 to 1979 and 1980 to 1999
comp = compare.periods(c("1960-01-01", "1979-12-31"), 
c("1980-01-01", "1999-12-31"), sycamore, plot=FALSE)

circ.s(comp)
}


\seealso{
\code{\link{compare.periods}}
}
