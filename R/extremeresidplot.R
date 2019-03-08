# special residual plots

residplot.extreme<-function(x,text=FALSE, data=FALSE) {
flows<-fourierAnalysis(x)

# get one sigma and two sigma events
s.hf<-sigmaHighFlows(x, resid.column=10)
s.lf<-sigmaLowFlows(x, resid.column=10)

# get annual extreme residuals
#syc.resids<-flows$signal$resid
resid.extremes<-annualExtremes(flows$signal, data.col=10,
                               moving.avg=FALSE)

max<-max(range(s.lf$onesigma.events$resid.sig)[2],
         range(s.hf$onesigma.events$resid.sig)[2])

min<-min(range(s.lf$onesigma.events$resid.sig)[2],
         range(s.hf$onesigma.events$resid.sig)[2])



#####################
## Plot 2 ###########
#####################

# plot of annual extreme residuals
# x axis is year, y axis is absolute residual

if(is.null(x$name)){
plot(resid.extremes$annual.max$year, abs(resid.extremes$annual.max$resid.sig), 
     pch=0, cex=1,
     xlab="Year", ylab="Residual event magnitude", 
     ylim=c(0, 1.15*max), main="Annual extreme residuals")
}
if(!is.null(x$name)){
  plot(resid.extremes$annual.max$year, abs(resid.extremes$annual.max$resid.sig), 
       pch=0, cex=1,
       xlab="Year", ylab="Residual event magnitude", 
       ylim=c(0, 1.15*max), main=paste("Annual extreme residuals:\n",x$name))
}
abline(h=2*s.lf$sigma.lfb, lty=2, col="lightskyblue3", lwd=2.5)

points(resid.extremes$annual.min$year, abs(resid.extremes$annual.min$resid.sig), 
       pch=1)
abline(h=2*s.hf$sigma.hfb, lty=3, col="seagreen", lwd=2.5)

if(text==TRUE) {
# the next lines add the year as text to the plot
text(s.hf$twosigma.events$year, s.hf$twosigma.events$resid.sig, labels=
       s.hf$twosigma.events$year, pos=2, cex=.7)

# the low flows are so concentrated that it looks bad for this data
text(s.lf$twosigma.events$year, abs(s.lf$twosigma.events$resid.sig), labels=
       s.lf$twosigma.events$year, pos=2, cex=.7)
}
if(data==TRUE){
return(resid.extremes)
}
}

