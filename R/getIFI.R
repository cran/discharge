# .....
# IFI
# .....
#' Inter-Flood Interval (IFI)
#' 
#' Compute Inter-Flood Interval (IFI) from the given residual values.
#' 
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param lowflow.start A vector giving start index of low-flow window in each year
#' @param lowflow.end A vector giving end index of low-flow window in each year
#' @param unique.years A vector or year values corresponding to the \code{highflow.start} 
#'                     and \code{highflow.end} values.
#' @param for.year (optional) Calculate IFI values only for the given year in this argument.
#'                 If argument is omitted, IFI values for all years are calculated. 
#' @return Data frame containing two columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{IFI} \tab Second column, represents IFI values
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
#' # signal parts
#' x.sp = getSignalParts(x.bl$pred2, candmin = c(40:125), candmax = c(190:330),
#'                       years = x.streamflow$data$year, 
#'                       months = x.streamflow$data$month,
#'                       jdays = x.streamflow$data$jday)
#' 
#' # IFI
#' ifi = getIFI(x.bl$resid.sig, x.streamflow$data$year, x.sp$LF.window.start, 
#'              x.sp$LF.window.end, x.sp$year)
#' 
#' @export
getIFI = function(resid, years, lowflow.start, lowflow.end, unique.years, for.year = NULL) {
  # validate inputs
  assert.for.year(for.year)
  assert.equal.length(resid, years)
  assert.equal.length(unique.years, lowflow.start, lowflow.end)
  
  if (is.null(for.year)) {
    lf.start = lowflow.start
    lf.end = lowflow.end
    window.years = unique.years
  } else {
    i = which(unique.years == for.year)
    lf.start = lowflow.start[i]
    lf.end = lowflow.end[i]
    window.years = unique.years[i]
  }
  
  ifi.year = c()
  ifi.value = c()
  for (i in seq_along(window.years)) {
    iyear = window.years[i]
    lf.start.thisyear = lf.start[i]
    lf.end.thisyear = lf.end[i]
    if (is.na(lf.start.thisyear)) {
      ifi.value = c(ifi.value, NA)
      ifi.year = c(ifi.year, iyear)
    } else {
      resid.window = resid[lf.start.thisyear:lf.end.thisyear]
      resid.window.rle = rle(resid.window > 0)
      resid.neg.lengths = resid.window.rle$lengths[which(resid.window.rle$values == FALSE)]
      ifi.value = c(ifi.value, max(0, resid.neg.lengths))
      ifi.year = c(ifi.year, iyear)
    }
  }
  
  ifi.data = data.frame(ifi.year, ifi.value)
  colnames(ifi.data) = c("year", "IFI")
  return(ifi.data)
}

