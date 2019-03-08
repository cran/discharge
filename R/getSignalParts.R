# .............
# Signal Parts
# .............
#' Signal parts
#' 
#' This function computes high flow and low flow window of seasonal signal, and the peak max
#' and peak min values.
#'
#' @param seas.sig Seasonal signal as generated from DFFT methods
#' @param candmin numeric vector of possible ordinal days in which the predicted signal is lowest. 
#'                This range need not be narrow, but a string of consecutive days should not include more than one local minimum. 
#'                Used for calculating the high- and low-flow windows.
#' @param candmax numeric vector of possible ordinal days in which the predicted signal is highest. 
#'                This range need not be narrow, but a string of consecutive days should not include more than one local maximum.
#' @param years A vector of years corrosponding to the seasonal signal values
#' @param months A vector of months corrosponding to the seasonal signal values
#' @param jdays A vector of julian days corrosponding to the seasonal signal values
#' @param for.year (optional) Calculate signal parts only for the given year in this argument.
#'                 If argument is omitted, all years are considered. 
#' @return Data frame containing following columns.
#' \tabular{ll}{
#' \code{year} \tab represents year \cr
#' \code{max.peak.index.all} \tab represents index value within the entire vector \cr
#' \code{max.peak.value} \tab represents value of max peak \cr
#' \code{highwind.start.index.all} \tab start index of high flow window within the entire vector \cr
#' \code{highwind.end.index.all} \tab end index of high flow window within the entire vector \cr
#' \code{lowwind.start.index.all} \tab start index of low flow window within the entire vector \cr
#' \code{lowwind.end.index.all} \tab end index of low flow window within the entire vector \cr
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
#' @export
getSignalParts = function(seas.sig, candmin, candmax, years, months, jdays, for.year = NULL) {
  # validate inputs
  assert.for.year(for.year)
  assert.equal.length(seas.sig, years, months, jdays)
  assert.numeric.vector(candmin)
  assert.numeric.vector(candmax)
  
  # Extract residual values and columns from the given series
  if (is.null(for.year)) {
    indices.years = c(1 : length(years))
    sig.years = years
  } else {
    indices.years = which(years == for.year)
    sig.years = years[indices.years]
  }
  
  if (!is.null(for.year) && ((for.year == years[1]) || (for.year == years[length(years)]))) {
    stop("'for.year' cannot be first or last year of the series")
  }
  
  unique.years = unique(sig.years)
  first.year = unique.years[1]
  last.year = unique.years[length(unique.years)]
  if (is.null(for.year)) {
    # exclude first and last years, as we analyse 20-month periods
    unique.years = tail(head(unique.years, -1), -1)
  }
  
  sigparts.year = c(first.year)
  sigparts.peak.index = c(NA)
  sigparts.peak.value = c(NA)
  sigparts.hfstart = c(NA)
  sigparts.hfend = c(NA)
  sigparts.lfstart = c(NA)
  sigparts.lfend = c(NA)
  for (iyear in unique.years) {
    # indices for this year
    indices.months.12 = which(years == iyear)
    # indices including last 4 months of past year and first 4 months of next year
    indices.months.20 = c(which((years == (iyear-1)) & (months > 8)), indices.months.12, which((years == (iyear+1)) & (months < 5)))
    
    # search for extrema of signals using candmin/candmax
    # minima
    indices.minblock=which(!is.na(match(jdays[indices.months.20], candmin))) # overlapping indices with candmin
    minblock.rle=rle(c(-10,indices.minblock)-c(indices.minblock,-10))
    count.minblock=length(which(minblock.rle$lengths>1)) # number of local minima to consider
    start.minblock=cumsum(minblock.rle$lengths)[which(minblock.rle$lengths>1)-1] # start of each block
    index.min=rep(NA,count.minblock)
    start.minblock=c(start.minblock,length(indices.minblock))
    
    for(j in 1:count.minblock){
      temp=indices.minblock[start.minblock[j]:(start.minblock[j+1]-1)]
      index.min[j]=temp[which.min(seas.sig[indices.months.20][temp])]
    }
    
    # maxima
    indices.maxblock=which(!is.na(match(jdays[indices.months.20], candmax))) # overlapping indices with candmax
    maxblock.rle=rle(c(-10,indices.maxblock)-c(indices.maxblock,-10))
    count.maxblock=length(which(maxblock.rle$lengths>1)) # number of local maxima to consider
    start.maxblock=cumsum(maxblock.rle$lengths)[which(maxblock.rle$lengths>1)-1] # start of each block
    
    index.max=rep(NA,count.maxblock)
    start.maxblock=c(start.maxblock,length(indices.maxblock))
    
    for(j in 1:count.maxblock){
      temp=indices.maxblock[start.maxblock[j]:(start.maxblock[j+1]-1)]
      index.max[j]=temp[which.max(seas.sig[indices.months.20][temp])]
    }
    
    # high and low flow windows
    if(length(index.max)>1) {
      lowflow.window.index=indices.months.20[index.max[1:2]]
    }
    
    if(length(index.min)>1) {
      highflow.window.index=indices.months.20[index.min[1:2]]
    }
    
    highflow.window=highflow.window.index[1]:highflow.window.index[2]
    lowflow.window=lowflow.window.index[1]:lowflow.window.index[2]
    
    ref.point.index = lowflow.window[which(years[lowflow.window.index]==iyear)][1] # first max peak in the specified year
    ref.point.value = seas.sig[ref.point.index]
    
    # append to vectors
    sigparts.year = c(sigparts.year, iyear)
    sigparts.peak.index = c(sigparts.peak.index, ref.point.index)
    sigparts.peak.value = c(sigparts.peak.value, ref.point.value)
    sigparts.hfstart = c(sigparts.hfstart, highflow.window.index[1])
    sigparts.hfend = c(sigparts.hfend, highflow.window.index[2])
    sigparts.lfstart = c(sigparts.lfstart, lowflow.window.index[1])
    sigparts.lfend = c(sigparts.lfend, lowflow.window.index[2])
  }
  
  # append NA values for the last year
  sigparts.year = c(sigparts.year, last.year)
  sigparts.peak.index = c(sigparts.peak.index, NA)
  sigparts.peak.value = c(sigparts.peak.value, NA)
  sigparts.hfstart = c(sigparts.hfstart, NA)
  sigparts.hfend = c(sigparts.hfend, NA)
  sigparts.lfstart = c(sigparts.lfstart, NA)
  sigparts.lfend = c(sigparts.lfend, NA)
  
  # prepare data frame
  sigparts.data = data.frame(sigparts.year, sigparts.peak.index, sigparts.peak.value, sigparts.hfstart, sigparts.hfend, sigparts.lfstart, sigparts.lfend)
  colnames(sigparts.data) = c("year", "peak.index", "peak.value", "HF.window.start", "HF.window.end", "LF.window.start", "LF.window.end")
  return(sigparts.data)
}

