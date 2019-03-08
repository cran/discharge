# .....
# HSAM
# .....
#' High Spectral Anomaly Mangitude (HSAM)
#' 
#' Compute High Spectral Anomaly Magnitude (HSAM) from the given residual values each year
#'
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param for.year (optional) Calculate HSAM values only for the given year in this argument.
#'                 If argument is omitted, HSAM values for all years are calculated.
#' @return Data frame containing four columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{HSAM} \tab Second column, represents HSAM values \cr
#' \code{index.year} \tab Third column, representing index of HSAM value in that \code{year} \cr
#' \code{index.all} \tab Fourth column, representing index of HSAM value in the input \code{resid} \cr
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
#' # HSAM
#' hsam = getHSAM(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' @export
getHSAM = function(resid, years, for.year = NULL) {
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
  hsam.year = c()
  hsam.value = c()
  hsam.index.year = c()
  hsam.index.all = c()
  for (iyear in unique.years) {
    indices.thisyear = which(resid.years == iyear)
    resid.values.thisyear = resid.values[indices.thisyear]
    index = which.max(resid.values.thisyear)
    hsam.year = c(hsam.year, iyear)
    hsam.value = c(hsam.value, resid.values.thisyear[index])
    hsam.index.year = c(hsam.index.year, index)
    hsam.index.all = c(hsam.index.all, indices.years[indices.thisyear[index]])
  }
  
  hsam.data = data.frame(hsam.year, hsam.value, hsam.index.year, hsam.index.all)
  colnames(hsam.data) = c("year", "HSAM", "Index.year", "Index.all")
  
  return(hsam.data)
}

