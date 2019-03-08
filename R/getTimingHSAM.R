# ............
# Timing HSAM
# ............
#' Time of occurence of High Spectral Anomaly Magnitude (HSAM)
#' 
#' Compute the number of days separating HSAM and reference point for each year.
#' 
#' @param index.hsam A scalar/vector of index of HSAM values in given year/years
#' @param index.ref A scalar/vector of index of reference point in given year/years
#' @param years A vector of years corresponding to HSAM and ref values.
#'               This argument can be NULL if the HSAM and ref values are scalars.
#' @param for.year (optional) Calculate timing (HSAM) only for the given year in this argument.
#'                 If argument is omitted, timing (HSAM) values for all years are calculated.
#' @return Scalar timing HSAM value if the inputs are scalars, or a Data frame containing two Columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{timing.hsam} \tab Second column, represents hsam timing values
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
#' # get HSAM values
#' hsam = getHSAM(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' # timing HSAM
#' thsam = getTimingHSAM(hsam$Index.all, x.sp$peak.index, x.sp$year)
#' 
#' @export
getTimingHSAM = function(index.hsam, index.ref, years, for.year = NULL) {
  # validate inputs
  assert.numeric.vector(index.hsam)
  assert.numeric.vector(index.ref)
  assert.numeric.vector(years)
  assert.equal.length(index.hsam, index.ref, years)
  assert.for.year(for.year)
  
  if (is.null(for.year)) {
    timing.hsam = abs(index.hsam - index.ref)
    timing.data = data.frame(years, timing.hsam)
  } else {
    indices.years = which(years == for.year)
    timing.hsam = abs(index.hsam[indices.years] - index.ref[indices.years])
    timing.data = data.frame(years[indices.years], timing.hsam)
  }
  colnames(timing.data) = c("year", "timing.hsam")
  return(timing.data)
}

