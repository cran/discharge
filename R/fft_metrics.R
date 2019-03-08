# .............
# DFFT Metrics
# .............
#' Discrete Fourier Transform Metrics
#' 
#' This is a wrapper function to calculate all the DFFT metrics for the given input signal
#' 
#' @param data A matrix with dates in the first column and discharge values in the second column.
#' Dates should be of the format "YYYY-MM-DD"
#' @param candmin numeric vector of possible ordinal days in which the predicted signal is lowest. This range need not be narrow, but a string of consecutive days should not include more than only local minimum. Used for calculating the high- and low-flow windows
#' @param candmax numeric vector of possible ordinal days in which the predicted signal is highest. This range need not be narrow, but a string of consecutive days should not include more than only local maximum.
#' @param river.name A character vector listing the river name.
#' @param baseline.signal If \code{NULL}, this function calculates baseline.signal using fourierAnalysis over the entire input series. The baseline signal can also be explicitly calculated and passed in as parameter. Check function \code{prepareBaseline()}
#' 
#' @return A list containing 2 data frames:
#' \tabular{ll}{
#' \code{high.level.metrics} \tab Data frame containing NAA and FPExt values for each year in the given series\cr
#' \code{naa.shape.components} \tab Data frame containing HSAM, LSAM, Transition time, HSAF, LSAF, timing of HSAM, timing of LSAM, IFI, IDI
#' }
#' 
#' @examples 
#' # load sample data
#' data("sycamore")
#' x = sycamore
#' 
#' # get streamflow object for the sample data
#' x.streamflow = asStreamflow(x)
#' 
#' # prepare baseline signal 
#' x.bl = prepareBaseline(x.streamflow)
#' 
#' # fetch the DFFT metrics for this sample data
#' # "candmax" chosen because preliminary analysis (e.g. with fourierAnalysis 
#' #           output) shows the signal is highest sometime between
#' #           day 190 and day 330
#' # "candmin" can be estimated analogously.
#' x.fftmetrics = fft_metrics(x, river.name = "Sycamore", candmin = c(40:125), 
#'                            candmax = c(190:330), baseline.signal = x.bl)
#' 
#' @export
fft_metrics = function(data, candmin, candmax, river.name = "", baseline.signal = NULL) {
  x.flows=asStreamflow(data, river.name=river.name)
  x.baseline = NULL
  x.resid = NULL
  x.years = x.flows$data$year
  x.months = x.flows$data$month
  x.jdays = x.flows$data$jday
  if (is.null(baseline.signal)) {
    x.sf = fourierAnalysis(x.flows)
    x.baseline = x.sf$signal$pred2
    x.resid = x.sf$signal$resid.sig
  } else {
    x.baseline = baseline.signal$pred2
    x.resid = baseline.signal$resid.sig
  }
  
  naa = getNAA(x.resid, x.years)
  fpext = getFPExt(x.resid, x.years)

  sp = getSignalParts(x.baseline, candmin=candmin, candmax=candmax, years=x.years, months=x.months, jdays=x.jdays)
  
  hsam = getHSAM(x.resid, x.years)
  lsam = getLSAM(x.resid, x.years)
  tt = getTransitionTime(hsam$Index.all, lsam$Index.all, hsam$year)
  hsaf = getHSAF(x.resid, x.years)
  lsaf = getLSAF(x.resid, x.years)
  thsam = getTimingHSAM(hsam$Index.all, sp$peak.index, sp$year)
  tlsam = getTimingLSAM(lsam$Index.all, sp$peak.index, sp$year)
  ifi = getIFI(x.resid, x.years, sp$LF.window.start, sp$LF.window.end, sp$year)
  idi = getIDI(x.resid, x.years, sp$HF.window.start, sp$HF.window.end, sp$year)
  
  high.level.metrics = cbind.data.frame(naa$year, naa$NAA, fpext$FPExt)
  colnames(high.level.metrics) <- c("year", "NAA", "FPExt")
  naa.shape.components = cbind.data.frame(naa$year, 
                                          hsam$HSAM, 
                                          lsam$LSAM, 
                                          tt$transition.time,
                                          hsaf$HSAF,
                                          lsaf$LSAF,
                                          thsam$timing.hsam,
                                          tlsam$timing.lsam,
                                          ifi$IFI,
                                          idi$IDI)
  colnames(naa.shape.components) <- c("year", "HSAM", "LSAM", "Transition.Time", "HSAM", "LSAF", "Timing.HSAM", "Timing.LSAM", "IFI", "IDI")
  
  return(list(high.level.metrics = high.level.metrics, naa.shape.components = naa.shape.components))
}

# data(etowah)
# fft.metrics = fft_metrics(etowah, candmin = c(190:330), candmax = c(40:125))
