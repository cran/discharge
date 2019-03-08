# .....
# IDI
# .....
#' Inter-Draught Interval (IDI)
#' 
#' Compute Inter-Draught Interval (IDI) from the given residual values.
#' 
#' @param resid A vector of residual values generated with respect to the baseline signal
#' @param years A vector of years corrosponding to the residual values
#' @param highflow.start A vector giving start index of high-flow window in each year
#' @param highflow.end A vector giving end index of high-flow window in each year
#' @param unique.years A vector or year values corresponding to the \code{highflow.start} 
#'                     and \code{highflow.end} values.
#' @param for.year (optional) Calculate IDI values only for the given year in this argument.
#'                 If argument is omitted, IDI values for all years are calculated. 
#' @return Data frame containing two columns:
#' \tabular{ll}{
#' \code{year} \tab First column, represents year \cr
#' \code{IDI} \tab Second column, represents IDI values
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
#' # IDI
#' idi = getIDI(x.bl$resid.sig, x.streamflow$data$year, x.sp$HF.window.start, 
#'              x.sp$HF.window.end, x.sp$year)
#' 
#' @export
getIDI = function(resid, years, highflow.start, highflow.end, unique.years, for.year = NULL) {
  # validate inputs
  assert.for.year(for.year)
  assert.equal.length(resid, years)
  assert.equal.length(unique.years, highflow.start, highflow.end)
  
  if (is.null(for.year)) {
    hf.start = highflow.start
    hf.end = highflow.end
    window.years = unique.years
  } else {
    i = which(unique.years == for.year)
    hf.start = highflow.start[i]
    hf.end = highflow.end[i]
    window.years = unique.years[i]
  }
  
  idi.year = c()
  idi.value = c()
  for (i in seq_along(window.years)) {
    iyear = window.years[i]
    hf.start.thisyear = hf.start[i]
    hf.end.thisyear = hf.end[i]
    if (is.na(hf.start.thisyear)) {
      idi.value = c(idi.value, NA)
      idi.year = c(idi.year, iyear)
    } else {
      resid.window = resid[hf.start.thisyear:hf.end.thisyear]
      resid.window.rle = rle(resid.window > 0)
      resid.pos.lengths = resid.window.rle$lengths[which(resid.window.rle$values == TRUE)]
      idi.value = c(idi.value, max(0, resid.pos.lengths))
      idi.year = c(idi.year, iyear)
    }
  }
  
  idi.data = data.frame(idi.year, idi.value)
  colnames(idi.data) = c("year", "IDI")
  return(idi.data)
}

