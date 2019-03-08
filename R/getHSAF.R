# .....
# HSAF
# .....
#' High Spectral Anomaly Frequency (HSAF)
#' 
#' Compute High Spectral Anomaly Frequency (HSAF) from the given residual values.
#' 
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param for.year (optional) Calculate HSAF values only for the given year in this argument.
#'                 If argument is omitted, HSAF values for all years are calculated. 
#' @return Data frame containing two Columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{HSAF} \tab Second column, represents HSAF values
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
#' # HSAF
#' hsaf = getHSAF(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' @export
getHSAF = function(resid, years, for.year = NULL) {
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
  hsaf.year = c()
  hsaf.value = c()
  for (iyear in unique.years) {
    # extract residual values
    resid.values.thisyear = resid.values[which(resid.years == iyear)]
    hsaf.rle = rle(resid.values.thisyear > 0)
    hsaf.value = c(hsaf.value, length(which(hsaf.rle$values == TRUE)))
    hsaf.year = c(hsaf.year, iyear)
  }
  
  hsaf.data = data.frame(hsaf.year, hsaf.value)
  colnames(hsaf.data) = c("year", "HSAF")
  return(hsaf.data)
}

