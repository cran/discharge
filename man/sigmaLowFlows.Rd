\name{sigmaLowFlows}
\alias{sigmaLowFlows}
\title{ Estimate catastrophic flow variability }
\description{
Calculates catastrophic variability for low flow events.  Negative residuals from the seasonal signal are used to calculate \eqn{\sigma.lf} , the standard deviation of low-flow events.
}
\usage{
sigmaLowFlows(x, resid.column)

}
\arguments{
\item{x}{ An object of class \code{data.frame} or \code{streamflow}.  If a \code{data.frame} is used, one column should contain residuals. }
\item{resid.column}{ Optional numeric specifiying which column contains residuals.  Required if \code{x} is a data frame.}

}
\value{
An object of class \code{list} with items
\item{n.droughts}{ Number of independent events with negative residuals. }
\item{sigma.lfa}{ Estimated sigma using the y-intercept. }
\item{sigma.lfb}{ Estimated sigma using the slope (\eqn{\sigma.lf}). }
\item{drought.line}{ Matrix containing fitted, observed, and residual values from regression of log counts on bin midpoints.  }
\item{onesigma.events}{ matrix containing information for all events below \eqn{\sigma.lf} (as calculated using the slope).  Columns will contain the same data as the output from the \code{\link{independentEvents}} function.  }
\item{twosigma.events}{ matrix containing information for all events below \eqn{2\sigma.lf}.  Columns will contain the same data as the output from the \code{\link{independentEvents}} function.  }
}

\examples{
# load data
data(sycamore)

# get streamflow object
sf = asStreamflow(sycamore)

# estimate catastrophic low flow variability
sigmaLowFlows(sf)
}

\seealso{
  \code{\link{independentEvents}}
  \code{\link{sigmaHighFlows}}
}

