\name{compare.periods}
\alias{compare.periods}
\alias{plot.compflows}
\title{ Compare residual variability across time periods }
\description{
The function finds the seasonal signal in the first of two time periods within the data series.  This signal is used to calculate residual flows in the first time period and "historical residuals", the residuals when the second period is compared to the signal in the first period. 

The output of the function gives event information for all residual events greater than \eqn{\sigma} as calculated from the functions \code{sigmaHighFlows} and \code{sigmaLowFlows}.
}
\usage{
compare.periods(p1, p2, x, plot=T)
\method{plot}{compflows} (x, ...)


}
\arguments{
\item{p1}{ Character vector specifying start and end date of the first period.  "YYYY-MM-DD" format.}
\item{p2}{ Character vector specifying start and end date of the second period.  "YYYY-MM-DD" format. }
\item{x}{ Matrix with first column specifying dates and second column specifying raw discharge data. }
\item{plot}{Logical; defaults to TRUE. If TRUE, the seasonal signal for both periods will be plotted.}
\item{...}{Other parameters. }
}
\value{
Object of S3 ckass \code{compflows} with the following items: 
\item{sigma.low}{Estimate of \eqn{\sigma-lf} from period 1 }
\item{sigma.high}{Estimate of \eqn{\sigma-hf} from period 1}

\item{p1.levents}{ Matrix of low flow events for period 1 }
\item{p1.hevents}{ Matrix of high flow events for period 1 }
\item{p2.levents}{ Matrix of low flow events for period 2 }
\item{p2.hevents}{ Matrix of high flow events for period 2 }

}

\examples{
# load data
data("sycamore")

# compare for periods from 1960 to 1979 and 1980 to 1999
compare.periods(c("1960-01-01", "1979-12-31"), 
c("1980-01-01", "1999-12-31"), sycamore)
}

\seealso{
\code{\link{sigmaHighFlows}}
\code{\link{sigmaLowFlows}}

}
