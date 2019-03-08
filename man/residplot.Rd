\name{residplot.extreme}
\alias{residplot.extreme}
\title{ Plot annual extreme residuals }
\description{
Creates a plot with the maximum annual low- and high residuals for each year in the series.}
\usage{
residplot.extreme(x, text=FALSE, data=FALSE)

}
\arguments{
\item{x}{ Object of class \code{streamflow}.}
\item{text}{ Logical.  If true, points corresponding to flows greater than \eqn{2\sigma} are labeled on the plot.}
\item{data}{Logical. If true, the extreme residuals are returned in the output.}
}
\value{
Plot with year on the x-axis and the maximum residual magnitude for that year on the y-axis.

If \code{data=TRUE}, output includes a \code{list} with the following components:

\item{annual.max}{ Matrix with data corresponding to the maximum residual flow for each year in series.}
\item{annual.min}{ Matrix with data corresponding to the minimum residual flow for each year in series.}

}

\examples{
\donttest{
# load data
data(sycamore)

# plot
residplot.extreme(asStreamflow(sycamore))
}
}


\seealso{
\code{\link{sigmaHighFlows}}
\code{\link{sigmaLowFlows}}

}
