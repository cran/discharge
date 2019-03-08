#........................
# Median finder (helper)
#........................
#' Find median value
#'
#' Find median value from the filtered part of the input vector.
#' This function is used as a helper to the bootstrapping routine in prepare baseline
#' 
#' @param data complete input vector
#' @param indices logical vector to filter indices from \code{data}
#' @return median value of the filtered vector
findMed = function(data, indices) {
  return(median(data[indices], na.rm=TRUE))
}


#..........
# baseline
#..........
#' Build baseline signal
#'
#' Runs fourier analysis on the input signal, to build baseline signal.
#' 
#' @param x streamflow object, as output from the \code{asStreamflow()} function
#' @param year.start Start of the year for estimating baseline, or \code{NULL} to interpret this from input data
#' @param year.end End of the year for estimating baseline, or \code{NULL} to interpret this from input data
#' @param window.20 If \code{TRUE}, baseline is constructed using windowing (20 year windows) and bootstrapping. 
#'                  If \code{FALSE}, baseline is constructed for a single run between start and end year.
#' @return \code{ssignal} object containing the baseline signal
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
#' bl.singlerun.all = prepareBaseline(x.streamflow)
#' 
#' # baseline for singlerun between the given start and end years
#' bl.singlerun.filtered = prepareBaseline(x.streamflow, year.start = 1961, 
#'                                          year.end = 2000)
#' 
#' # baseline with windowinng and bootstrapping on all years in the input signal
#' bl.windowed.all = prepareBaseline(x.streamflow, window.20 = TRUE)
#' 
#' # baseline with windowing and bootstrapping on given start year 
#' #  with end year inferred from singal
#' bl.windowed.filtered = prepareBaseline(x.streamflow, year.start = 1961, 
#'                                        window.20 = TRUE)
#' 
#' @export
prepareBaseline = function(x, year.start = NULL, year.end = NULL, window.20 = FALSE) {
  
  if (class(x) != "streamflow") {
    stop("x must be a stream flow object")
  }
  
  # Extract rows from x, corresponding to the range of years specified
  if (is.null(year.start)) year.start = x$data$year[1]
  if (is.null(year.end)) year.end = x$data$year[length(x$data$year)]
  data = x$data
  x.filterobj = x
  filtered.indices = which((data$year >= year.start) & (data$year <= year.end))
  x.filtered = x$data[filtered.indices,]
  x.filterobj$data = x.filtered
  
  # work on length of data
  n.tot = dim(x.filtered)[1]
  even <- n.tot%%2==0
  if (!even) {x.filtered <- x.filtered[-n.tot,]}
  n.tot = dim(x.filtered)[1]
  whole.years = floor(n.tot/365)
  x.filtered = x.filtered[1:(whole.years*365),]
  
  # decision branch... windowed or complete analysis
  is.seasonal = TRUE
  pred2.366days = rep(0, 366)
  if (window.20) {
    # prepare baseline signal for each 20 year consecutive window
    window.length = 20
    window.end = year.end - window.length + 1
    total.windows = year.end - year.start - window.length + 2
    signal.matrix = matrix(rep(0.0, total.windows*366),nrow = 366, ncol = total.windows)
    window.obj = x.filterobj
    for (window.index in year.start:(year.start + total.windows - 1)) {
      # extract window indices and update object data
      window.indices = which(x.filtered$year >= window.index & x.filtered$year <= (window.index+window.length-1))
      window.obj$data = x.filtered[window.indices,]
      # detrend the window data
      index = 1:length(window.obj$data$ldis.corrupt)
      fit.index = lm(window.obj$data$ldis.corrupt ~ index)
      window.obj$data$ldis.corrupt = fit.index$residuals
      
      # perform fourier analysis with stationary = true (data already detrended)
      window.sf = fourierAnalysis(window.obj, stationary = T)
      
      # append detrended reconstructed signal to signals matrix (column wise)
      signal.matrix[,window.index - year.start + 1] = window.sf$signal$pred2[1:366]
    }
    
    # for each day, resample baseline datapoints and select median from the resampled data
    bootobj.list = apply(signal.matrix, 1, boot, findMed, 1)
    
    for (i in 1:366) {
      pred2.366days[i] = bootobj.list[[i]]$t0
    }

  } else {
    # detrend, and apply fourier analysis with stationary = T
    index = 1:length(x.filterobj$data$ldis.corrupt)
    fit.index = lm(x.filterobj$data$ldis.corrupt ~ index)
    x.filterobj$data$ldis.corrupt = fit.index$residuals
    
    sf = fourierAnalysis(x.filterobj)
    pred2.366days = sf$signal$pred2[1:366]
  }
  
  # Reconstruction: Add trend to the detrended signal
  # Extract trend for the entire signal
  index = 1:length(data$ldis.corrupt)
  fit.index = lm(data$ldis.corrupt ~ index)
  signal.trend = fit.index$fitted.values
  
  # repeat seasonal (detrended) signal for the entire range of years
  jday = data$jday
  signal.detrended = rep(0, length(index))
  for (i in 1:366) {
    signal.detrended[jday == i] = pred2.366days[i]
  }

  # add trend to the detrended sinusoidal baseline
  pred2 = signal.trend + signal.detrended

  # compute residuals
  resid.sig = data$ldis.corrupt - pred2
  
  # extend the data matrix
  signal = as.data.frame(cbind(data, resid.sig, pred2))
  
  return(signal)
}

# # source("C:/Users/Samarth/Documents/ASU/FutureH2O/discharge_0.1.1/discharge/R/testfunctions.R")
# data(sycamore)
# sycamore.flows=asStreamflow(sycamore, river.name="Sycamore")
# yr = 1980
# x.sf = fourierAnalysis(sycamore.flows)
# 
# # x.fftmetrics = fftmetrics(sycamore.flows,yr,candmin=c(190:330),candmax=c(40:125))
# 
# x.baseline = prepareBaseline(sycamore.flows, window.20 = FALSE)
# 
# ann.ind=c(which(x.baseline$year==(yr-1) & x.baseline$month>8),which(x.baseline$year==yr),which(x.baseline$year==(yr+1) & x.baseline$month<5))
# plot(x.baseline[1:10000,8],cex=.75,
#      ylab="log normalized discharge",xlab="ordinal day",xaxt="n")
# #days=x.baseline$signal[ann.ind,7]
# #d.seq=seq(1,595,45)+1
# #axis(1,at=seq(1,595,45)+1,as.character(days[d.seq]))
# points(x.baseline[1:10000,11],type="l",lwd=3,col="blue")
