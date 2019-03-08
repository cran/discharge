\name{sigmaHighFlows}
\alias{sigmaHighFlows}
\title{ Estimate catastrophic flow variability }
\description{
Calculates catastrophic variability for high flow events.  Positive residuals from the seasonal signal are used to calculate \eqn{\sigma.hf} , the standard deviation of high-flow events.
}
\usage{
sigmaHighFlows(x, resid.column)

}
\arguments{
\item{x}{ An object of class \code{data.frame} or \code{streamflow}.  If a \code{data.frame} is used, one column should contain residuals. }
\item{resid.column}{ Optional numeric specifiying which column contains residuals.  Required if \code{x} is a data frame.}

}
\value{
An object of class \code{list} with items
\item{n.floods}{ Number of independent events with positive residuals. }
\item{sigma.hfa}{ Estimated sigma using the y-intercept. }
\item{sigma.hfb}{ Estimated sigma using the slope (\eqn{\sigma.hf}). }
\item{flood.line}{ Matrix containing fitted, observed, and residual values from regression of log counts on bin midpoints.  }
\item{onesigma.events}{ matrix containing information for all events below \eqn{\sigma.hf} (as calculated using the slope).  Columns will contain the same data as the output from the \code{\link{independentEvents}} function.  }
\item{twosigma.events}{ matrix containing information for all events below \eqn{2\sigma.hf}.  Columns will contain the same data as the output from the \code{\link{independentEvents}} function.  }
}

\examples{
# load data
data(sycamore)

# get streamflow object
sf = asStreamflow(sycamore)

# estimate catastrophic high flow variability
sigmaHighFlows(sf)
}


\seealso{
  \code{\link{independentEvents}}
  \code{\link{sigmaLowFlows}}
}