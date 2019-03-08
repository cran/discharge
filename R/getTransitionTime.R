# ................
# Transition Time
# ................
#' Transition Time
#' 
#' Compute the number of days separating HSAM and LSAM for the given year/years.
#' 
#' @param index.hsam A scalar/vector of index of HSAM values in given year/years
#' @param index.lsam A scalar/vector of index of LSAM values in given year/years
#' @param years A vector of years corresponding to HSAM and LSAM values.
#'               This argument can be NULL if the HSAM and LSAM values are scalars.
#' @param for.year (optional) Calculate transition time only for the given year in this argument.
#'                 If argument is omitted, transition times for all years are calculated. 
#' @return Scalar transition time if the inputs are scalars, or a Data frame containing two Columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{transition.time} \tab Second column, represents transition times
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
#' # get HSAM and LSAM values
#' hsam = getHSAM(x.bl$resid.sig, x.streamflow$data$year)
#' lsam = getLSAM(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' # transition time
#' tt = getTransitionTime(hsam$Index.all, lsam$Index.all, hsam$year)
#' 
#' @export
getTransitionTime = function(index.hsam, index.lsam, years, for.year = NULL) {
  # validate inputs
  assert.numeric.vector(index.hsam)
  assert.numeric.vector(index.lsam)
  assert.equal.length(index.hsam, index.lsam, years)
  assert.for.year(for.year)
  
  if (is.null(for.year)) {
    transition_time = abs(index.hsam - index.lsam)
    trx.data = data.frame(years, transition_time)
  } else {
    indices.years = which(years == for.year)
    transition_time = abs(index.hsam[indices.years] - index.lsam[indices.years])
    trx.data = data.frame(years[indices.years], transition_time)
  }
  colnames(trx.data) = c("year", "transition.time")
  return(trx.data)
}

