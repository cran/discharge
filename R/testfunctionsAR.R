annualExtremes<-function(x, data.col=NULL, year.col=NULL, moving.avg=FALSE)
{
  
  
  
  if(class(x)=="streamflow") data<-x$data
  if (class(x)=="data.frame" | class(x)=="matrix") data<-x
  if(is.null(data)) stop("Data in unrecognized format; try dataframe or matrix.")
  n<-dim(data)[1]
  
  
  
  if (is.null(data.col)) {
    if (class(x)!="streamflow") stop ("Cannot process matrix or dataframe without data.col")
    
    if (moving.avg==FALSE)
    {
      data.vec<-data$ldis.corrupt }
    if (moving.avg==TRUE)
    {
      data.vec<-data$smoothed.7d
    }
    
  }
  
  if (!is.null(data.col) & !is.numeric(data.col)) {
    
    stop("data.col should be column number of data.")}
  if (is.numeric(data.col)) {
    data.vec<-data[,data.col]
    if(length(data.vec)!=n) warning("data.vec is different length than data matix")
    
  }
  
  if(is.null(data$year) & is.null(year.col)) stop ("Missing year.vec")
  if(is.null(data$year) & !is.null(year.col)) year.vec<-data[,year.col]
  if (!is.null(data$year)) year.vec<-data$year
  
  
  index<-as.integer(1:length(data.vec))
  year<-unique(year.vec)
  n.years<-length(year)
  
  
  max<-NULL
  min<-NULL
  
  for (i in 1:n.years)
  {
    
    data.thisyear<-which(data$year==year[i])
    
    if (length(which(is.na(data.vec[data.thisyear])==FALSE))<1)
    {
      warning ("all obs are NA in year ", year[i])
      
      break
    }
    if (length(data.thisyear)==1) warning("Only one obs in year ", year[i])
    
    max.datathisyear<-data.thisyear[which(data.vec[data.thisyear]==
                                            max(data.vec[data.thisyear], na.rm=T))]
    n.max<-length(max.datathisyear)
    max<-c(max,index[max.datathisyear])
    
    
    min.datathisyear<-data.thisyear[which(data.vec[data.thisyear]==
                                            min(data.vec[data.thisyear], na.rm=T))]
    n.min<-length(min.datathisyear)
    min<-c(min,index[min.datathisyear])
  }
  
  
  annual.max<-data[max,]
  
  annual.min<-data[min,]
  
  
  return(list(annual.max=annual.max, annual.min=annual.min))
}



annualnoise<-function(x)
{

  if (class(x)=="streamflow")
  {

    annual.max<-annualExtremes(x)$annual.max$ldis.corrupt
  }
  if (class(x)!="numeric" & class(x)!="streamflow") stop ("Argument x must be numeric vector 
                                   or streamflow object.")
  if(class(x)=="numeric"){
  annual.max<-x}
  n.max<-length(annual.max)
  power.of.two<-c(2,4,8,16,32,64,128,256, 512, 1024)
  M.ind<-max(which(power.of.two<=n.max))
  num.lags<-power.of.two[M.ind]
  
  auto.corr.list<-acf(annual.max, type="correlation", demean="TRUE",
                      lag.max=num.lags, main="ACF", plot=FALSE, na.action=na.pass)
   auto.corr<-auto.corr.list$acf[-1]
  
  #omit first obs because zero lag will always have corr of 1 #
  max.corr<-max(auto.corr)
  
  lag.max.corr<-which(auto.corr==max.corr)
  order<-lag.max.corr
  

  lci<-qnorm(.025)/sqrt(num.lags)
  uci<--lci
  nyquist<-num.lags/2
  
  significant<-which(abs(auto.corr)>uci)
  order<-significant[significant==lag.max.corr]
  
  if (length(order)>0)
  {
    auto.corr[significant]
    
    fit.ar<-ar(annual.max, aic=FALSE, order.max=order, method="yule-walker")
    phi.params<-fit.ar$ar 
    
    #transform to frequency domain
    ddt<-order
    freq<-c((0:nyquist)/(2*nyquist*ddt))
    ddd<-matrix(data=NA, ncol=order, nrow=length(freq))
    den.ps<-rep(NA, length(freq))
    for (j in 1:length(freq))
    {
      for (i in 1:order)
      {
        ddd[j,i]<-phi.params[i]*exp(-2*pi*i*sqrt(as.complex(-1))*freq[j])
      }
      den.ps[j]<-(abs(1-sum(ddd[j,])))^2
    }
    
    # check on this
    w2<-rep(var(auto.corr, na.rm=TRUE), length(freq))
    power.spec.acf<-w2/den.ps
                                                              
    log.ps<-log(power.spec.acf[2:length(power.spec.acf)])
    log.freq<-log(freq[2:length(freq)])
    vv<-as.data.frame(cbind(log.freq, log.ps))
    
    #determine noise color by regression of log power on log freq
    fit.vlog<-lm(log.ps~log.freq, data=vv)
    
    int<-coefficients(fit.vlog)[1]
    # noise color
    theta.a<-coefficients(fit.vlog)[2]
    theta.se<-coef(summary(fit.vlog))[2,]
    df.error<-fit.vlog$df.residual
    ci.up<-theta.a+theta.se*qt(.975, df=df.error)
    ci.low<-theta.a-theta.se*qt(.975, df=df.error)
    ci.theta.a<-cbind(ci.up, ci.low)
    
    if( (ci.low<=0 && ci.up >=0)) { #insignificant slope
      int<-mean(log.ps)
      theta.a<-0
      
      warning("note:order is not zero, but regression parameters not significant")
    } 
  } else { #loop for when no significant ACF
    #warning("note: order is zero")
    fit.ar<-NULL
    ddt<-1
    freq<-c((0:nyquist)/(2*nyquist*ddt))
    order<-0
    ppp<-(var(annual.max, na.rm=TRUE))*exp(-2*pi*sqrt(as.complex(-1)) )
    abs.ppp<-(abs(ppp))^2
    power.spec<-rep(abs.ppp, length(freq))
    fit.vlog<-lm(log(abs.ppp)~1)
    int<-coefficients(fit.vlog)[1]
    theta.a<-0
    log.freq<-log(freq)
    vv<-cbind(log.freq,log(power.spec))
  }
  
  rho.txm<-as.data.frame(cbind(1:length(auto.corr), round(auto.corr, digits=4)))
  names(rho.txm)<-c("lag", "autocorr")
  interval<-cbind(lci, uci)

  
  output<-list(auto.corr=rho.txm, lm.fit=fit.vlog, 
               interval=interval, log.log=vv, reg.stats=c(int, theta.a),order=order,  fit.ar=fit.ar)
  class(output)<-"annualnoise"
  return(output)
  
}

print.annualnoise<-function(x)
{
  cat("Autoregressive model of order", x$order, ". \n ")
  print(x$fit.ar)
  cat("\n Sample autocorrelations: \n \n")
  print(x$auto.corr)
  cat("\n \n")
  cat("Noise color:", x$reg.stats[2], "\n")
  if (!(x$order==0))
  {
   cat("Regression of log power spectrum on log frequency: \n")
  print(summary( x$lm.fit))
  }
  
}

summary.annualnoise<-function(object, ...)
{
  if (!is.null(object$fit.ar)) {
    cat("Autoregressive model of order", object$order, ". \n \n")
  print(object$fit.ar) }
  if (is.null(object$fit.ar)) {cat("No significant autocorrelation.  \n")}
  cat("\n \n")
  cat("Noise color:", object$reg.stats[2], "\n")

}


lp3Events<-function(x)
{

  
  if (class(x)=="streamflow") {
  annual.max<-annualExtremes(x, moving.avg=FALSE)$annual.max$ldis.corrupt
  annual.min<-annualExtremes(x,moving.avg=FALSE)$annual.min$ldis.corrupt
  
  annual.7dmax<-annualExtremes(x, moving.avg=TRUE)$annual.max$smoothed.7d
  annual.7dmin<-annualExtremes(x,moving.avg=TRUE)$annual.min$smoothed.7d
  }
  if(class(x)!="matrix" & class(x)!="streamflow") stop ("Input must be matrix or streamflow object.")
  if (class(x)=="matrix") {
  annual.max<-annualExtremes(x[,2], moving.avg=FALSE)$annual.max$ldis.corrupt
  annual.min<-annualExtremes(x[,1], moving.avg=FALSE)$annual.max$ldis.corrupt
  annual.7dmax<-annualExtremes(x[,2], moving.avg=TRUE)$annual.max$smoothed.7d
  annual.7dmin<-annualExtremes(x[,1],moving.avg=TRUE)$annual.min$smoothed.7d
  }

  annual.max<-as.numeric(annual.max)
  annual.min<-as.numeric(annual.min)

 
  # stores log pearson III parameters
  params<-pelpe3(samlmu(annual.max))
  
  #ten year and two year events
  Q10<-as.numeric(quape3(.90, params))
  BFD<-as.numeric(quape3(.50, params))
  
  params.low<-pelpe3(samlmu(annual.7dmin))
  #ten year and two year events
  d7.low10<-as.numeric(quape3(.10, params.low))
  d7.low2<-as.numeric(quape3(.50, params.low))
  
  event.def<-list(Q2=BFD, Q10=Q10, L2=d7.low2, L10=d7.low10)
  
  return(event.def)
}



