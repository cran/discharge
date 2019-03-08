fourierAnalysis<-function(x, stationary=F)
{ 
  if (!is.null(x$na.info)) {
  n.blocks<-dim(x$na.info)[1]
  warning (paste(n.blocks, "blocks of more than 10 nas.  Results may be inaccurate."))
  
  }
  
  data<-x$data
  n.tot<-dim(data)[1]
  even<-n.tot%%2==0
  if (!even) {data<-data[-n.tot,]}
  n.tot<-dim(data)[1]
  whole.years<-floor(n.tot/365)
  data<-data[1:(whole.years*365),]
  jdayfirst<-data[1,7]
  n.tot<-dim(data)[1]
  
  n.trend<-n.tot
  na.obs<-which(is.na(data$ldis.corrupt))
  not.na<-which(!is.na(data$ldis.corrupt))
  
  #detrend data by regressing on index
  index<-1:n.tot
  fit.index <- lm(data$ldis.corrupt ~ index)
  predicted<-predict(fit.index, newdata=as.data.frame(index))
  detrended.data<-rep(NA, n.tot)
  if (stationary==F) {
  detrended.data[not.na]<-(data$ldis.corrupt-predicted)[not.na]
  detrended.data [na.obs]<- 0
  int.roi <- coefficients(fit.index)[1]
  slope.roi <- coefficients(fit.index)[2]
  }
  if (stationary==T) {
  detrended.data<-data$ldis.corrupt
  detrended.data[na.obs]<-predicted[na.obs]
  int.roi<-mean(detrended.data, na.rm=T)
  slope.roi<-0 
  warning(paste("Stationary model used. Estimated trend has coefficient",
  round(coefficients(fit.index)[2], 6), "p value", 
  round(summary(fit.index)$coefficients[2,4], digits=7)))
  }

  
  n.nas<-length(which(is.na(detrended.data)))
  
  dt<-1/365
  
  full.freq<-(seq(0, n.trend/2, by=1) ) / (n.trend*dt)
  n.fr<-length(full.freq)
  index.flowfreq<-which(full.freq==1.0 | full.freq==1.5 | full.freq==2 |
                          full.freq==3 | full.freq==4 | full.freq==(1/3))
  
  full.fr<-rep(0, n.tot)
  full.fr[1:n.fr]<-full.freq
  
  filler.l<-length(index.flowfreq)
  if (length(filler.l)==0)
    stop("No exact match to frequencies")
  
  #fourier transform, phase and amplitude
  fft<-fft(detrended.data)
  amp.spec<-Mod(fft)/(length(detrended.data)/2)
  power.spec<-amp.spec^2
  phase<-Arg(fft)
  phase<-phase-(((jdayfirst-1)*2*pi*full.fr)/365)
  
  #modifies to give phase in range -pi to pi
  for (i in 1:length(phase))
  {while (phase[i]< (-1*pi)){
  phase[i]<-phase[i]+2*pi}}
  
  N<-n.trend
  m<-floor(N/5)
  
  # Hanning window 
  all.hann.amps<-c(.5*(power.spec[2]+power.spec[3]),.25*(power.spec[2:(n.trend-2)]
                                                         +2*power.spec[3:(n.trend-1)]+power.spec[4:n.trend]), .5*(power.spec[n.trend-1]+
                                                                                                                    power.spec[n.trend]))
  
  hann.amps<-c(.5*(power.spec[2]+power.spec[3]),.25*(power.spec[2:(m-2)]
                                                     +2*power.spec[3:(m-1)]+power.spec[4:m]), .5*(power.spec[m-1]+
                                                                                                    power.spec[m]))
  
  nu<-2.67*(N/m) #degrees of freedom for the chi-square
  
  seqData<-cbind(hann.amps[index.flowfreq], seq(1, filler.l, by=1))
  seqData2<-seqData[order(seqData[,1], decreasing=TRUE),] 
  seqData3<-cbind(seqData2, seq(1, filler.l, by=1))
  seqData4<-seqData3[order(seqData3[,2]),]
  
  lci<-rep(NA, filler.l)
  uci<-rep(NA, filler.l)
  
  for (i in 1:filler.l)
  {
    no.peaks.tested<-seqData4[i,3]
    index.temp<-index.flowfreq[i]
    
    ci.lower<-nu*hann.amps/qchisq((1-.025)/no.peaks.tested, df=nu)
    ci.upper<-nu*hann.amps/qchisq(.025/no.peaks.tested, df=nu)
    
    l.ci.index<-ci.lower[index.temp-1]
    u.ci.index.leftn<-ci.upper[index.temp-4]
    u.ci.index.rightn<-ci.upper[index.temp+2]
    u.ci.index.avg<-max(u.ci.index.leftn, u.ci.index.rightn)
    med.u.ci<-median(ci.upper)
    
    lci[i]<-l.ci.index
    uci[i]<-u.ci.index.avg
    
  }
  
  sig.amp.index<-which(lci>uci)
  
  ## log log scale for regression
  new.amps<-c(999, all.hann.amps)
  log.power<-log(new.amps)[2:((length(new.amps)/2)+1)] #take log and shorten
  log.freq<-log(full.freq[2:length(full.freq)])
  # quality check
  if (length(log.power)!=length(log.freq))
    stop("error: log.power and new.freq are different lengths")
  
  if (length(sig.amp.index)>0)
  {
    seasonal<-1
    significant.flowfreq<-(index.flowfreq[sig.amp.index])
    n.significant<-length(significant.flowfreq)
    target.amps<-amp.spec[significant.flowfreq]
    target.freqs<-full.freq[significant.flowfreq]
    target.phase<-phase[significant.flowfreq]
    idds<-seq(1, length(significant.flowfreq))
    terms<-cbind(target.freqs, target.amps, target.phase, idds)
  }
  
  if(length(sig.amp.index)==0)
  {
    seasonal<-0
    terms<-0
  }
  #regression on log-log scale to get theta #
  #remove significant frequencies #
  if (seasonal==1)
  {
    new.log.power<-log.power[-significant.flowfreq] #this is output at the end
    new.log.freq<-log.freq[-significant.flowfreq] #this is output at the end
  } else {
    new.log.power<-log.power
    new.log.freq<-log.freq
  }
  
  x1<-rep(1, length(new.log.freq))
  xx<-as.data.frame(cbind(x1, new.log.freq, new.log.power))
  
  lm.powerspec<-lm(new.log.power~new.log.freq, data=xx)
  
  rms<-findRMS(new.log.power, terms, seasonal)
  
  #get stats from regression fit
  stats<-summary(lm.powerspec)
  r.squared<-stats$r.squared
  theta<--1*coefficients(lm.powerspec)[2]
  level<-coefficients(lm.powerspec)[1]
  fitted.power<-fitted(lm.powerspec)
  
  if (seasonal==1)
  {
    stats<-list(r.squared=r.squared, theta=theta, level=level, terms=terms,
                slope.roi=slope.roi, int.roi=int.roi) #.roi means stats from regression on index
    logpower<-cbind(log.power, log.freq)
    
    ## log of power / freq with significant signal removed.
    noise.data<-as.data.frame(cbind(new.log.power, new.log.freq, fitted.power))
    
    # call predictSignalResid to construct signal
    signal<-predictSignalResid(terms, data, int.roi, slope.roi)
    obj<-list("terms"=terms, "signal"=signal, "detrend.fit"=fit.index,
    "logps.regression"=lm.powerspec, seasonal=seasonal, 
    log.powerfreq=logpower, rms=rms,river.name=x$name)
  } else {
    
    terms<-list(paste("no significant signal")   )
    pred2<-rep(int.roi, length(data$ldis.corrupt)) #intercept is used for signal
    resid.sig<-data$ldis.corrupt-pred2
    signal<-as.data.frame(cbind(data, resid.sig, pred2))
    logpower<-cbind(log.power, log.freq)
    obj<-list("terms"=terms, "signal"=signal, "detrend.fit"=fit.index,
    "logps.regression"=lm.powerspec,  seasonal=seasonal, 
    log.powerfreq=logpower, rms=rms, river.name=x$name)

  }
  
  class(obj)<-"ssignal"
  
  
  
  return(obj)
}

