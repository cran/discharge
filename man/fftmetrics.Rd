\name{fftmetrics}
\alias{fftmetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of spectral anomolies for a specified year.}
\description{
This function takes as input a streamFlow object and a year, and outputs the timing and magnitude of noteworthy spectral anomolies.
}
\usage{
fftmetrics(x,year,candmin,candmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{streamflow object, as output from the asStreamflow function.
}
\item{year}{integer; data from this year will be analyzed.  The plotted hydrograph will include data from August of the previous year to April of the following year.}

\item{candmin}{numeric vector of possible ordinal days in which the predicted signal is lowest.  This range need not be narrow, but a string of consecutive days should not include more than only local minimum. Used for calculating the high- and low-flow windows.}
\item{candmax}{numeric vector of possible ordinal days in which the predicted signal is highest.  This range need not be narrow, but a string of consecutive days should not include more than only local maximum.}
}
\value{
A list object with the following components;
\item{sam}{dataframe where row 1 corresponds to the largest low-residual event in the year, and row 2 corresponds to the largest high-residual.  The dataframe contains each event's date, ordinal day, magnitude of residual and signal, index in the original data, and timing (in days) relative to the reference point.}
\item{ref.point}{date used as a reference point for the timing of max and min events.  If not given by the user, it is the first local maximum of the signal, within the year specified.  }
\item{events}{dataframe with rows corresponding to start dates of the longest low-flow and high-flow events.  The data frame contains signal and residual data for the start date and the lenght of the run.}
\item{auc}{dataframe with three values: 'net.auc', 'rel.auc.low', and rel.auc.high'. 'net.auc' is the sum of positive residuals in the high flow window divided by the sum of negative residuals in the low flow window. 'rel.auc.low' is postive residuals divided by the sum of negative residuals }
\item{noise.color}{numeric, "theta.d" as calculated as in the fourierAnalysis function.}
}

\examples{
data(etowah)
etowah.flows=asStreamflow(etowah, river.name="Etowah")
# "candmax" chosen because preliminary analysis (e.g., with fourierAnalysis output) 
# shows the signal is highest sometime between day 40 and day 125.
# "candmin" can be estimated analogously.
fftmetrics(etowah.flows,2002,candmin=c(190:330),candmax=c(40:125)) }