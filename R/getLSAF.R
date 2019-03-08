# .....
# LSAF
# .....
#' Low Spectral Anomaly Frequency (LSAF)
#' 
#' Compute Low Spectral Anomaly Frequency (LSAF) from the given residual values.
#' 
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param for.year (optional) Calculate LSAF values only for the given year in this argument.
#'                 If argument is omitted, LSAF values for all years are calculated. 
#' @return Data frame containing two Columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{LSAF} \tab Second column, represents LSAF values
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
#' # LSAF
#' lsaf = getLSAF(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' @export
getLSAF = function(resid, years, for.year = NULL) {
  # validate inputs
  assert.for.year(for.year)
  assert.equal.length(resid, years)
  
  if (is.null(for.year)) {
    resid.values = resid
    resid.years = years
  } else {
    indices.years = which(years == for.year)
    resid.values = resid[indices.years]
    resid.years = years[indices.years]
  }
  
  unique.years = unique(resid.years)
  lsaf.year = c()
  lsaf.value = c()
  for (iyear in unique.years) {
    # extract residual values
    resid.values.thisyear = resid.values[which(resid.years == iyear)]
    lsaf.rle = rle(resid.values.thisyear > 0)
    lsaf.value = c(lsaf.value, length(which(lsaf.rle$values == FALSE)))
    lsaf.year = c(lsaf.year, iyear)
  }
  
  lsaf.data = data.frame(lsaf.year, lsaf.value)
  colnames(lsaf.data) = c("year", "LSAF")
  return(lsaf.data)
}

