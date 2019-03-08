sigmaLowFlows<-function(x, resid.column)
{
  
  if (class(x)=="streamflow")  {
    fft<-fourierAnalysis(x)
    data<-fft$signal
    resid.column<-10
  }
  
  
  if (class(x)=="matrix" | class(x)=="data.frame" ) {
    if( !is.null(resid.column)  )
      data<-x
    if (is.null(resid.column)) stop("specify column of residuals")
  }
  # call independent events function with residuals as data vec
  r2<-independentEvents(cutoff.val=0, data, resid.column, below.cutoff=TRUE)
  
  # duration of more than 1 day
  r2<-r2[which(r2[,3]>1),]
  
  minresids.event<-r2[,4]
  n.droughts<-length(minresids.event)
  
  
  range<-range(minresids.event, na.rm=TRUE)
  s<-seq(range[1], 0, by=(abs(range[1])/5))
  h<-hist(minresids.event, breaks=s, plot=FALSE)
  
  log.counts<-log(h$counts+1)
  mids.sq<-h$mid^2
  ones<-rep(1, length(log.counts))
  
  # fit half normal model
  x3<-as.data.frame(cbind(log.counts, ones, mids.sq))
  fit<-lm(log.counts~mids.sq, data=x3)
  beta.0<-coefficients(fit)[1]
  beta.1<-coefficients(fit)[2]
  
  beta.0.mult<-coef(summary(fit))[1,]*qt(.025, df=4)
  beta.1.mult<-coef(summary(fit))[2,]*qt(.025, df=4)
  
  beta.0.int<-c(beta.0-beta.0.mult, beta.0, beta.0+beta.0.mult)
  beta.1.int<-c(beta.1-beta.1.mult, beta.1, beta.1+beta.1.mult)
  
  x1<-mids.sq
  predicted<-fitted(fit)
  observed<-log.counts
  residuals<-resid(fit)
  drought.line<-as.data.frame(cbind(x1, predicted, observed, residuals))
  
  p.value<-coef(summary(fit))[2,4]
  
  sigma.lf.a<-(1/(sqrt(2*pi))*n.droughts/exp(beta.0)) #times 2?
  sigma.lf.b<-sqrt(1/(-beta.1*2))
  
  # get events less than 1 sigma
  onesigma.events<-independentEvents(cutoff.val=-sigma.lf.b, data, data.column=resid.column,
                                     below.cutoff=TRUE)                                                        
  # call independent events to find events less than 2sigma
  twosigma.events<-independentEvents(cutoff.val=-2*sigma.lf.b, data, data.column=resid.column,
                                     below.cutoff=TRUE) 
  
  return(list(n.droughts=n.droughts, sigma.lfa=sigma.lf.a, sigma.lfb=sigma.lf.b, 
              drought.line=drought.line, onesigma.events=onesigma.events,
              twosigma.events=twosigma.events))
  
} #end else loop



sigmaHighFlows<-function(x, resid.column)
{
	

  if (class(x)=="streamflow")  {
   fft<-fourierAnalysis(x)
   data<-fft$signal
   resid.column<-10
   }


   if (class(x)=="matrix" | class(x)=="data.frame" ) {
   if( !is.null(resid.column)  )
   data<-x
   if (is.null(resid.column)) stop("specify column of residuals")
   }

r2<-independentEvents(cutoff.val=0.00000001, data, resid.column, below.cutoff=FALSE)

r2<-r2[which(r2[,3]>1),]

maxresids.event<-r2[,4]
n.floods<-length(maxresids.event)

range<-range(maxresids.event, na.rm=TRUE)
s<-seq(0, range[2]+.001, by=(abs(range[2])/5))
h<-hist(maxresids.event, breaks=s, plot=FALSE)


log.counts<-log(h$counts+1)
mids.sq<-h$mid^2
ones<-rep(1, length(log.counts))

x2<-as.data.frame(cbind(log.counts, ones, mids.sq))
fit<-lm(log.counts~mids.sq, data=x2)
beta.0<-coefficients(fit)[1]
beta.1<-coefficients(fit)[2]

beta.0.mult<-coef(summary(fit))[1,]*qt(.025, df=4)
beta.1.mult<-coef(summary(fit))[2,]*qt(.025, df=4)

beta.0.int<-c(beta.0-beta.0.mult, beta.0, beta.0+beta.0.mult)
beta.1.int<-c(beta.1-beta.1.mult, beta.1, beta.1+beta.1.mult)

x1<-(mids.sq) 
predicted<-fitted(fit)
observed<-log.counts
residuals<-resid(fit)
flood.line<-as.data.frame(cbind(x1, predicted, observed, residuals))

p.value<-coef(summary(fit))[2,4]

sigma.hf.a<-(1/(sqrt(2*pi))*n.floods/exp(beta.0)) #times 2?
sigma.hf.b<-sqrt(1/(-beta.1*2))

onesigma.events<-independentEvents(cutoff.val=sigma.hf.b, data, resid.column,
                                   below.cutoff=FALSE) 

twosigma.events<-independentEvents(cutoff.val=2*sigma.hf.b, data, resid.column,
                                   below.cutoff=FALSE) 

return(list(n.floods=n.floods, sigma.hfa=sigma.hf.a, sigma.hfb=sigma.hf.b, 
 flood.line=flood.line, p.value=p.value, onesigma.events=onesigma.events, 
                        twosigma.events=twosigma.events))

} #end function loop


