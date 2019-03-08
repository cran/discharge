\name{sigmaeventsplot}
\alias{sigmaeventsplot}
\title{ Plot events by day of the year }
\description{
Creates a .  This is similar to the plots produced in the "compare.periods" function, but only displays the data for a single time period.

}
\usage{
sigmaeventsplot(x)


}
\arguments{
\item{x}{ Object of class "streamflow".}
}
\value{
A "ggplot2" plot depicting frequency of events greater than sigma, organized circularly by ordinal day of the year.
}

\seealso{
\code{\link{sigmaHighFlows}}
\code{\link{sigmaLowFlows}}
\code{\link{compare.periods}}

}
