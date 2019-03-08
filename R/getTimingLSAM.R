# ............
# Timing HSAM
# ............
#' Time of occurence of Low Spectral Anomaly Magnitude (LSAM)
#' 
#' Compute the number of days separating LSAM and reference point for each year.
#'
#' @param index.lsam A scalar/vector of index of LSAM values in given year/years
#' @param index.ref A scalar/vector of index of reference point in given year/years
#' @param years (optional) A vector of years corresponding to LSAM and ref values.
#'               This argument can be NULL if the LSAM and ref values are scalars.
#' @param for.year (optional) Calculate timing (LSAM) only for the given year in this argument.
#'                 If argument is omitted, timing (LSAM) values for all years are calculated.
#' @return Scalar timing LSAM value if the inputs are scalars, or a Data frame containing two Columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{timing.lsam} \tab Second column, represents lsam timing values
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
#' # get signal parts
#' x.sp = getSignalParts(x.bl$pred2, candmin = c(40:125), candmax = c(190:330),
#'                       years = x.streamflow$data$year, 
#'                       months = x.streamflow$data$month,
#'                       jdays = x.streamflow$data$jday)
#' 
#' # get LSAM values
#' lsam = getLSAM(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' # timing LSAM
#' tlsam = getTimingLSAM(lsam$Index.all, x.sp$peak.index, x.sp$year)
#' 
#' @export
getTimingLSAM = function(index.lsam, index.ref, years = NULL, for.year = NULL) {
  # validate inputs
  assert.numeric.vector(index.lsam)
  assert.numeric.vector(index.ref)
  assert.numeric.vector(years)
  assert.equal.length(index.lsam, index.ref, years)
  assert.for.year(for.year)
  
  if (is.null(for.year)) {
    timing.lsam = abs(index.lsam - index.ref)
    timing.data = data.frame(years, timing.lsam)
  } else {
    indices.years = which(years == for.year)
    timing.lsam = abs(index.lsam[indices.years] - index.ref[indices.years])
    timing.data = data.frame(years[indices.years], timing.lsam)
  }
  colnames(timing.data) = c("year", "timing.lsam")
  return(timing.data)
}


