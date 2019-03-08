# ....
# NAA
# ....
#' Net Annual Anomaly (NAA)
#' 
#' Calculate Net Annual Anomaly (NAA) from the given residual values.
#' 
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param for.year (optional) Calculate NAA values only for the given year in this argument.
#'                 If argument is omitted, NAA values for all years are calculated. 
#' @return Data frame containing two columns: 
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{NAA} \tab Second column, represents NAA values
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
#' # NAA
#' naa = getNAA(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' @export
getNAA = function(resid, years, for.year = NULL) {
  # check inputs
  assert.for.year(for.year)
  assert.equal.length(resid, years)
  
  # Extract residual values and columns from the given series
  if (is.null(for.year)) {
    resid.values = resid
    resid.years = years
  } else {
    indices.years = which(years == for.year)
    resid.values = resid[indices.years]
    resid.years = years[indices.years]
  }
  
  unique.years = unique(resid.years)
  naa.year = c()
  naa.value = c()
  for (iyear in unique.years) {
    # extract residual values
    resid.values.thisyear = resid.values[which(resid.years == iyear)]
    naa.year = c(naa.year, iyear)
    naa.value = c(naa.value, sum(resid.values.thisyear, na.rm = TRUE))
  }
  
  naa.data = data.frame(naa.year, naa.value)
  colnames(naa.data) = c("year", "NAA")
  
  return(naa.data)
}
