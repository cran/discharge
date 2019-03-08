# .....
# LSAM
# .....
#' Low Spectral Anomaly Mangitude (LSAM)
#' 
#' Compute Low Spectral Anomaly Magnitude (LSAM) from the given residual values each year
#'
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param for.year (optional) Calculate LSAM values only for the given year in this argument.
#'                 If argument is omitted, LSAM values for all years are calculated.
#' @return Data frame containing four columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{LSAM} \tab Second column, represents LSAM values \cr
#' \code{index.year} \tab Third column, representing index of LSAM value in that \code{year} \cr
#' \code{index.all} \tab Fourth column, representing index of LSAM value in the input \code{resid} \cr
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
#' # LSAM
#' lsam = getLSAM(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' @export
getLSAM = function(resid, years, for.year = NULL) {
  # validate inputs
  assert.for.year(for.year)
  assert.equal.length(resid, years)
  
  # Extract residual values and columns from the given series
  if (is.null(for.year)) {
    indices.years = c(1 : length(years))
    resid.values = resid
    resid.years = years
  } else {
    indices.years = which(years == for.year)
    resid.values = resid[indices.years]
    resid.years = years[indices.years]
  }
  
  unique.years = unique(resid.years)
  lsam.year = c()
  lsam.value = c()
  lsam.index.year = c()
  lsam.index.all = c()
  for (iyear in unique.years) {
    indices.thisyear = which(resid.years == iyear)
    resid.values.thisyear = resid.values[indices.thisyear]
    index = which.min(resid.values.thisyear)
    lsam.year = c(lsam.year, iyear)
    lsam.value = c(lsam.value, resid.values.thisyear[index])
    lsam.index.year = c(lsam.index.year, index)
    lsam.index.all = c(lsam.index.all, indices.years[indices.thisyear[index]])
  }
  
  lsam.data = data.frame(lsam.year, lsam.value, lsam.index.year, lsam.index.all)
  colnames(lsam.data) = c("year", "LSAM", "Index.year", "Index.all")
  
  return(lsam.data)
}


