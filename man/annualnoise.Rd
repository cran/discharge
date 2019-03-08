\name{annualnoise}
\alias{annualnoise}
\alias{summary.annualnoise}
\title{ Annual noise color }
\description{ 
 The autocovariance function is estimated for the annual maxima in the series.  An autoregressive model of the order of the highest significant lag is fit, using the Yule-Walker method to estimate the parameters.  The function is transformed into the frequency domain, yielding an estimate \code{theta.a} of the annual noise color.
}
\usage{
annualnoise(x)
}
\arguments{
\item{x}{ A numeric vector of annual extremes. A \code{streamflow} object may also be used.  If input is \code{streamflow} the function uses annual maximum discharge. }
}
\value{
An object of S3 class \code{annualnoise} with the following attributes:
\item{auto.corr}{ Sample autocorrelation. }
\item{lm.fit}{ \code{lm} object from regression of log power spectrum on log frequency. }
\item{interval}{ Upper and lower bounds of a 95\% acceptance region when \eqn{\rho=0}.  }
\item{log.log}{ Matrix with log frequency and log power spectrum. }
\item{reg.stats}{Slope and intercept of regression of log power spectrum on log frequency, where slope is the annual noise color (\code{theta.a}).}
\item{order}{ Indicates order of fitted AR model. }
\item{fit.ar}{Object of class \code{ar} summarizing the fitted AR model.}

}
\details{
To determine the order of the AR model, the ACF is calculated at all lags less than or equal to the highest power of 2 less than the length of the series.  The order of the AR model is the lag with the highest significantly non-zero autocorrelation.
}
\examples{
data(sycamore)
sycamore.flows<-asStreamflow(sycamore,river.name="Sycamore Creek")
syc.ar<-annualnoise(sycamore.flows)
}
