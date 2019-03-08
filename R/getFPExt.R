# ......
# FPExt
# ......
#' Flood Pulse Extent (FPExt)
#' 
#' Calculate the Flood Pulse Extent (FPExt) from the given residual values.
#' 
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param for.year (optional) Calculate FPExt values only for the given year in this argument.
#'                 If argument is omitted, NAA values for all years are calculated. 
#' @return Data frame containing two columns: 
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{FPExt} \tab Second column, represents FPExt values
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
#' # FPExt
#' fpext = getFPExt(x.bl$resid.sig, x.streamflow$data$year)
#' 
#' @export
getFPExt = function(resid, years, for.year = NULL) {
  # validate inputs
  assert.for.year(for.year)
  assert.equal.length(resid, years)
  
  # Extract residual values and columns from the given series
  resid.values = NULL
  resid.years = NULL
  if (is.null(for.year)) {
    resid.values = resid
    resid.years = years
  } else {
    indices.years = which(years == for.year)
    resid.values = resid[indices.years]
    resid.years = years[indices.years]
  }
  
  unique.years = unique(resid.years)
  fpext.year = c()
  fpext.value = c()
  for (iyear in unique.years) {
    # extract residual values
    resid.values.thisyear = resid.values[which(resid.years == iyear)]
    mean.thisyear = mean(resid.values.thisyear)
    fpext.year = c(fpext.year, iyear)
    fpext.value = c(fpext.value, 
                  sum(resid.values.thisyear[which(resid.values.thisyear >= mean.thisyear)]
                                 , na.rm = TRUE))
  }
  
  fpext.data = data.frame(fpext.year, fpext.value)
  colnames(fpext.data) = c("year", "FPExt")
  
  return(fpext.data)
}
