


subsetData<-function(range.start, range.end, data.matrix, max.nas=10)
{


Range<-c(range.start, range.end)
range.date<-as.Date(Range, "%Y-%m-%d")
range<-as.numeric(as.Date(Range, "%Y-%m-%d"))
range.year<-as.numeric(format(range.date, format="%Y"))
n.years<-(range.year[2]-range.year[1])
max.year<-max(as.numeric(data.matrix$year))-1

sub.data<-data.matrix[data.matrix$date.num>=range[1] &
data.matrix$date.num<=range[2],]
l.span<-max.year-n.years

finite.logical<-is.finite(sub.data[,2])
n.nas<-length(finite.logical[finite.logical==FALSE])
if (n.nas>0)
{
if (n.nas>max.nas)
stop("Set has ", n.nas, " NA values", sep=" ")
else
warning("data have ", n.nas, " NA values", sep=" ")
}

return(list("sub.data"=sub.data, "l.span"=l.span, n.nas=n.nas))

}




predictSignalResid<-function(dterms, data.matrix, int, slope)
{


fund.freq<-dterms[,1]
fund.amp<-dterms[,2]
fund.phase<-dterms[,3]

ldis.corrupt<-data.matrix[,8]
jday<-data.matrix[,7]

d.p<-seq(1:length(data.matrix[,1]))/365
index.all<-seq(1:length(jday))
l<-seq(1:366)
l.terms<-length(dterms[,1])
seq<-seq(1, l.terms, by=1)
piece<-matrix(data=NA, nrow=length(d.p), ncol=l.terms)
sig.sep<-matrix(data=NA, nrow=length(l), ncol=l.terms)
sig.com<-rep(NA, length(l))

for (i in 1:l.terms)
{
piece[,i]<-fund.amp[i]*cos(2*pi*d.p*fund.freq[i]+fund.phase[i])
}

for (i in 1:length(piece[,1]))
{
sig.com[i]<-sum(piece[i,])+int
}

xs<-seq(1:length(jday))*slope  
pred2<-rep(NA, length(jday))
for (i in 1:length(jday))
{
pred2[i]<-xs[i]+sig.com[jday[i]]
}

resid.sig<-ldis.corrupt-pred2
data.matrix<-as.data.frame(cbind(data.matrix, resid.sig, pred2))

return(data.matrix)
}


## s3 methods for ssignal
summary.ssignal<-function(object, ...)
{
cat("Signal info", "\n")
print(object$terms)
cat("\n","Noise color:", -1*coefficients(object$logps.regression)[2], "\n")
cat("\n","Average discharge:", coefficients(object$detrend.fit)[1], "\n")
cat("\n","Signal-noise ratio:", object$rms$snr, "\n")
}



plot.ssignal<-function(x, plot.type="hydrograph", ...)
{
if (plot.type=="hydrograph")
{
  river.name<-x$river.name
plot(x$signal$jday, x$signal$ldis.corrupt,
cex=.35, main=paste("Seasonal signal:", river.name),
xlab="ordinal day", ylab="log discharge")
jdayone<-x$signal$jday[1]
dayshift<-365-jdayone
pred.ind<-(dayshift+1):(dayshift+366)
lines(1:366, x$signal$pred2[pred.ind], col="red", lwd=4)
}

if (plot.type=="autocorr")
{

d<-x$signal
resids<-d$resid.sig
corr.coefs<-acf(resids, lag.max=365, na.action=na.pass)$acf
plot(1:366, corr.coefs, type="l", main="Daily autocorrelation", xlab="lags",
ylab="correlation coefficient", ylim=c(-1,1))
xx<-c(1:365, 365:1)
yy<-c(corr.coefs[1:365], rep(-1, 365))
polygon(xx, yy, col="black")

xloc<-366
yloc<-.95
text(x=xloc, y=yloc, paste(river.name), pos=2)

}


}

print.ssignal<-function(x, ...)
{

print(list(x$terms, detrend.lm=x$detrend.fit,
 power.spec.lm=x$logps.regression, seasonal=as.logical(x$seasonal)), ...)
print(x$rms, ...)

}


## function to find RMS
###################################################
##### findRMS #####################################
###################################################

findRMS<-function(new.log.power, dterms, seasonal)
  
{
  

  amps.long<-sqrt(exp(new.log.power)) 
  
  
  if (seasonal==0)
  {
    a<-sum(amps.long^2)
    rms.noise<-sqrt(a/length(amps.long))
    t<-list(rms.noise=rms.noise, "no significant signal")
  } else {
    if (dim(dterms)[1]==1) {
      sig.amps<-dterms[2]} else {
      sig.amps<-dterms[,2]}
    rms.signals<-sqrt(mean(sig.amps^2))
    
    a<-sum(amps.long^2)
    b<-sum(sig.amps^2)
    
    rms.noise<-sqrt( a/length(amps.long) )
    snr<-20*log(rms.signals/rms.noise)
    
    t<-list(rms.signal=rms.signals, rms.noise=rms.noise,
            ratio=rms.signals/rms.noise, snr=snr)
  }
  return(t)
}

independentEvents<-function(cutoff.val, data, data.column, below.cutoff=FALSE)
{
  
  if( is.null(data.column)) stop ("Arg data.column missing.")
  
  data.vec<-data[,data.column]
  index.vec<-1:(dim(data)[1])
  
  mult.for.nas<--1
  if (below.cutoff==TRUE)
    mult.for.nas<-1
  
  seq<-seq(1,length(data.vec), by=1)
  na.index<-which(is.na(data.vec==TRUE))
  #warning commented for now.  maybe work on suppressing it later
  #if (length(na.index)>0)
  #  warning("NAs in data.vec. considered non event days.")
  data.vec[na.index]<-mult.for.nas*999  #replaces NA values with 999 or -999

  
  if (below.cutoff==TRUE)
  {
    dev<-data.vec<=cutoff.val
    extreme.direction<-"min"
    sort.direction<-"FALSE"
  } else {
    dev<-data.vec>=cutoff.val
    extreme.direction<-"max"
    sort.direction<-"TRUE"
    mult.for.nas<-1
  }
  
  
  run.lengths<-rle(dev)
  trues<-which(run.lengths$value==TRUE)
  length.of.trues<-as.numeric(run.lengths$length[trues])
  first.event<-run.lengths$values[1] 
  event.duration<-length.of.trues
  
  #end function if there are no true events
  if (length(trues)==0) warning (paste("No events meeting criterion.",extreme.direction, cutoff.val))
  
  #long loop for when there is at least one event
  if (length(trues)>0)
  {
    n.events<-length(trues)
    events.starts<-rep(NA, n.events)
    events.ends<-rep(NA, n.events)
    extreme.this.event<-rep(0, n.events)
    ind.extreme<-rep(NA, n.events)
    duplicates<-rep(NA, n.events)
    
    #loop for when first run is false, rle output needs to be handled differently
    if (!first.event)
    {
     for (i in 1:(length(trues)))
     {
       start<-sum(run.lengths$lengths[1:(2*i-1)])+1
       end<-start+length.of.trues[i]-1
       events.starts[i]<-start
       events.ends[i]<-end
       event.data<-data.vec[start:end]
       event.index<-index.vec[start:end]
       
       if (sort.direction=="TRUE")
       {
         extreme.thisevent<-which(event.data==max(event.data))
       } else {
         extreme.thisevent<-which(event.data==min(event.data))
         
       }
    
       #warning for multiple extremes#
       if (length(extreme.thisevent)>1)
         duplicates[i]<-1
       if (length(extreme.thisevent)==1)
         duplicates[i]<-0
       
       extreme.this.event[i]<-event.data[extreme.thisevent[1]]
       ind.extreme[i]<-event.index[extreme.thisevent[1]] #lowest index
     } 
     }  else { #begin loop for first run true
      

      events.starts[1]<-1
      events.ends[1]<-1+length.of.trues[1]-1
      
      if (events.ends[1]!=1) {
        event.data<-data.vec[1:events.ends[1]]
        event.index<-index.vec[1:events.ends[1]]
                                    } else {
        event.data<-data.vec[1]
        event.index<-index.vec[1]
                                    }
        if (sort.direction=="TRUE")
        {
          extreme.thisevent<-which(event.data==max(event.data))
        } else {
          extreme.thisevent<-which(event.data==min(event.data))
          
        } 
        
        if (length(extreme.thisevent)>1)
          duplicates[1]<-1
        if (length(extreme.thisevent)==1)
          duplicates[1]<-0
        
        extreme.this.event[1]<-event.data[extreme.thisevent[1]]
        ind.extreme[1]<-event.index[extreme.thisevent[1]] 
        
            
        for (i in 2:(length(trues)))
        {
          start<-sum(run.lengths$lengths[1:(2*i-2)])+1
          end<-sum(run.lengths$lengths[1:(2*i-1)])
          
          events.starts[i]<-start
          events.ends[i]<-end
          
          event.data<-data.vec[start:end]
          event.index<-index.vec[start:end]
          
          if (sort.direction=="TRUE")
          {
            extreme.thisevent<-which(event.data==max(event.data))
          } else {
            extreme.thisevent<-which(event.data==min(event.data))
            
          }
    
          if (length(extreme.thisevent)>1)
            duplicates[i]<-1
          if (length(extreme.thisevent)==1)
            duplicates[i]<-0
          
          extreme.this.event[i]<-event.data[extreme.thisevent[1]]
          ind.extreme[i]<-event.index[extreme.thisevent[1]] 
          
          
          
        } #end true-start loop
      }
    }      
    
  if (length(trues)>0) {
    matrix<-as.data.frame(cbind(events.starts, events.ends, event.duration,
                               extreme.this.event, ind.extreme, data[ind.extreme,],
                               duplicates)) }
    
    return(matrix)
    
  
}
