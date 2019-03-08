\name{fourierAnalysis}
\alias{fourierAnalysis}
\alias{plot.ssignal}
\alias{summary.ssignal}
\title{Extract seasonal signal from time series }
\description{
  The fast Fourier transform is used to extract the seasonal signal of a time series.  The significant frequencies are found from among periods of length 2-, 3-, 4-, 6-, 12-, and 18-months.  

The signal may be specified as stationary or non-stationary.  If a non-stationary fit is allowed, simple linear regression estimates the long term linear trend.   The seasonal signal is calculcated from the residuals.

Predicted flow (and corresponding residual) at each time point is calculated from seasonal signal and, if non-stationary, long term trend coefficient.

}
\usage{
fourierAnalysis(x, stationary=F)
\method{plot}{ssignal} (x, plot.type="hydrograph", ...)
}
\arguments{
\item{x}{ An object of class \code{streamflow}}
\item{stationary}{ Logical; defaults to FALSE.}
\item{plot.type}{ Indicates the type of plot to create.  The default "hydrograph" produces a plot of ordinary day and log normalized discharge, with the seasonal signal overlaid.  "auto.corr" produces a plot of daily autocorrelation as calculated from the residual flows.}
\item{...}{Other parameters. }
}
\value{
An object of class \code{ssignal} with items
\item{signal}{ Data matrix augmented to included predicted and residual values. }
\item{terms}{ Matrix containing amplitude, phase, and frequency of seasonal signal.  }
\item{detrend.fit}{ An \code{lm} object from regression of discharge on index of observation.  }
\item{logps.regression}{ An \code{lm} object from regression of log power spectrum on log frequency (where log frequencies have seasonal signal removed.) }
\item{rms}{ \code{list} containing RMS amplitude for noise, RMS amplitude for signal, and signal-to-noise ratio. }
}
\examples{
data(sycamore)
sycamore.flows<-asStreamflow(sycamore,river.name="Sycamore Creek")
syc.seas<-fourierAnalysis(sycamore.flows)
summary(syc.seas)
}

