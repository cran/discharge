#........................
# Baseline filter
#........................
#' Filter the baseline signal for a given time window
#'
#' @param bl baseline signal as returned from the function \code{prepareBaseline()}
#' @param filter.date.start start date of the filtering window
#' @param filter.date.end end date of the filtering window
#' @param date.format format of date specified in \code{filter.date.start} and \code{filter.date.end}
#' @return baseline signal filtered for the given date window
#' 
#' @examples 
#' # load sample data
#' data("sycamore")
#' x = sycamore
#' 
#' # get streamflow object for the sample data
#' x.streamflow = asStreamflow(x)
#' 
#' # baseline for single run for all the years in input signal
#' x.bl = prepareBaseline(x.streamflow)
#' 
#' # filter the baseline signal between years 1993 and 2000
#' x.bl.filtered = filterBaseline(x.bl, filter.date.start = "1993-01-01",
#'                                filter.date.end = "2000-12-31")
#' 
#' @export
filterBaseline = function(bl, filter.date.start, filter.date.end, date.format = "%Y-%m-%d") {
  dates = as.Date(bl$date, format=date.format)
  start.date = as.Date(filter.date.start, format=date.format)
  end.date = as.Date(filter.date.end, format=date.format)
  bl.filtered = bl[which(dates >= start.date & dates <= end.date),]
  return(bl.filtered)
}
