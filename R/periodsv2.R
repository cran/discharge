

compare.periods<-function(p1, p2, x, plot=T) {

 # validation
  
if (length(p1)!=2) stop ("p1 must vector of two dates")
if (length(p2)!=2) stop ("p2 must vector of two dates")

name.dummy<-NULL
# match dates

p1.data<-asStreamflow(x,name.dummy,start.date=p1[1], end.date=p1[2], max.na=5000)
p2.data<-asStreamflow(x,name.dummy,start.date=p2[1], end.date=p2[2],max.na=5000)

p1.flows<-fourierAnalysis(p1.data, stationary=T)
p2.flows<-fourierAnalysis(p2.data, stationary=T)

if( plot=="T") {
par(ask=T)
plot(p1.flows)
plot(p2.flows)
}
p1.flowdata<-p1.flows$signal
p2.flowdata<-p2.flows$signal


# want a matrix with jday and prediction
pred.1<-p1.flowdata$pred2[which(p1.flowdata$jday==1)]
if (length(unique(pred.1))!=1) stop ("non-stationary predictions")

start.day<-which(p1.flowdata$jday==1)[1]
first.366<-which(p1.flowdata$jday==366)[1] 
unique.predictions<-p1.flowdata$pred2[start.day:(start.day+364)]
unique.predictions<-c(unique.predictions, p1.flowdata$pred2[first.366])

#find historical residuals
firstday.p2<-p2.flowdata$jday[2]
n.p2<-dim(p2.flowdata)[1]
p2.pred<-rep(NA, n.p2)
for (j in 1:n.p2) {
  p2.pred[j]<-unique.predictions[p2.flowdata$jday[j]]
}
hist.resid<-p2.flowdata$ldis.corrupt-p2.pred
p2.flowdata<-cbind(p2.flowdata,hist.resid)

# find unique pred for period 2
# want a matrix with jday and prediction
pred.2<-p2.flowdata$pred2[which(p2.flowdata$jday==1)]
if (length(unique(pred.2))!=1) stop ("non-stationary predictions")

start.day2<-which(p2.flowdata$jday==1)[1]
first.366.2<-which(p2.flowdata$jday==366)[1] 
unique.predictions2<-p2.flowdata$pred2[start.day2:(start.day2+364)]
unique.predictions2<-c(unique.predictions2, p2.flowdata$pred2[first.366.2])

if( plot==T) {
  yrange<-c(min(c(unique.predictions,unique.predictions2))*.95,
            max(c(unique.predictions,unique.predictions2))*1.1)
plot(unique.predictions, type="l", main="Fitted signal: period 1 and period 2",
     ylab="Predicted", xlab="Ordinal day", ylim=yrange)
lines(unique.predictions2, lty=3)
legend("topright", legend=c("Period 1", "Period 2"), lty=c(1,3))
}
par(ask=F)

# find sigma from first period
sigma.lf.p1<-sigmaLowFlows(p1.flowdata, 10)
sigma.hf.p1<-sigmaHighFlows(p1.flowdata, 10)

slf.p1<-sigma.lf.p1$sigma.lfb
shf.p1<-sigma.hf.p1$sigma.hfb

p2.onesigmal.events<-independentEvents(-slf.p1,p2.flowdata, 12, below.cutoff=T)
p2.onesigmah.events<-independentEvents(shf.p1,p2.flowdata, 12, below.cutoff=F)

data.frames<-list(sigma.low=slf.p1, sigma.high=shf.p1,p1.levents=sigma.lf.p1$onesigma.events,
                  p1.hevents=sigma.hf.p1$onesigma.events, 
                  p2.levents=p2.onesigmal.events,
                  p2.hevents=p2.onesigmah.events)
class(data.frames)<-"compflows"

return(data.frames)
}

getbins<-function(range){
  
  rres<-range
  rres[1]<-rres[1]-.0001
  rres[2]<-rres[2]+.0001
  
  if (.75<(rres[2]-rres[1])) dig<-1 else dig<-2
  if(rres[1]>0){
  decimal<-decimalplaces(abs(signif(rres[1],1)-signif(rres[1],2)))}
  if(rres[1]<0){
    decimal<-decimalplaces(abs(signif(rres[2],1)-signif(rres[2],2)))}
  start<-(floor(rres[1]*10^decimal)/10^decimal)
  end<-ceiling(rres[2]*10^decimal)/10^decimal
  int<-(end-start+.01)/5
  if(rres[1]>0) {
    return(seq(start,start+5*int,int))   }
  if(rres[1]<0) {
    return(seq(start+5*int,start,-int))   }
}


plot.compflows<-function(x, ...){
  if(class(x)!="compflows") stop ("Input should be output from compare.periods.")
  .e<-environment()
  lf.1<-x[[3]]
  lf.1$jday[which(lf.1$jday==366)]<-365
  hf.1<-x[[4]]
  hf.1$jday[which(hf.1$jday==366)]<-365
  lf.2<-x[[5]]
  lf.2$jday[which(lf.2$jday==366)]<-365
  hf.2<-x[[6]]
  hf.2$jday[which(hf.2$jday==366)]<-365
  
  circ<-round(circ.s(x),3)
  lf.1.c<-circ[2,2:4]
  hf.1.c<-circ[1,2:4]
  hf.2.c<-circ[3,2:4]
  lf.2.c<-circ[4,2:4]
  
  ################################
  # these define breaks for x and where to put ticks on plot
  q.test<-seq(0,365,14)
  q.test[length(q.test)]<-366
  x.ticks<-c(0,31,59,90,120,151,181,212,243,273,304,334)
  label.v2<-c("J","F","M","A","M","J","J","A","S","O","N","D")
  
  #low flow p1
  rres<-c(min(range(lf.1$resid.sig)[1],range(lf.2$hist.resid)[1]),
          max(range(lf.1$resid.sig)[2],range(lf.2$hist.resid)[2]))
  
  cuts<-getbins(rres)
  residual.magnitude1<-factor(cut(lf.1$resid.sig,
                                  cuts),ordered=T)
  t1<-reorder(residual.magnitude1, -as.numeric(residual.magnitude1))
  residual.magnitude3<-cut(lf.2$hist.resid,
                           cuts)
  t3<-reorder(residual.magnitude3, -as.numeric(residual.magnitude3))
  # to set y range for low flow
  wtable1<-table(residual.magnitude1, cut(lf.1$jday,q.test))
  wtable3<-table(residual.magnitude3, cut(lf.2$jday,q.test))
  yrange<-range(c(as.numeric(wtable1),as.numeric(wtable3)))
  if(yrange[2]>=15){
    max.count<-ceiling(yrange[2]/5)
    ymaxl<-max.count*5 +5
  }
  if(yrange[2]<15){
    max.count<-ceiling(yrange[2]/3)
    ymaxl<-max.count*3 +3
  }
  if (ymaxl<=15) by.int<-3
  if (ymaxl<30 & ymaxl>15) by.int<-5
  if (ymaxl>=30) by.int<-10
  
  
  plot.1<-ggplot(lf.1, aes(x=lf.1$jday,fill=t1), environment = .e)+
    geom_histogram(breaks=q.test) +
    scale_y_continuous(expand=c(0,0),limits= c(0,ymaxl),breaks=seq(0,ymaxl,by.int))+
    scale_x_discrete(breaks= x.ticks,labels=label.v2)+
    coord_polar() + ggtitle("Period 1, Low residuals") +
    scale_fill_discrete(name=paste("Summary statistics:\nmu=", bquote(.(lf.1.c[1])), "\nrho=", bquote(.(lf.1.c[2])),
                                   "\nkappa=", bquote(.(lf.1.c[3])),
                                   
                                   "\n\n\n\nResidual magnitude"))+
    xlab("Ordinal day")
  
  suppressMessages(print(plot.1))
  
  
  
  
  #high flow p1
  rres<-c(min(range(hf.1$resid.sig)[1],range(hf.2$hist.resid)[1]),
          max(range(hf.1$resid.sig)[2],range(hf.2$hist.resid)[2]))
  cuts<-getbins(rres)
  residual.magnitude2<-factor(cut(hf.1$resid.sig,
                                  cuts), ordered=T)
  residual.magnitude4<-factor(cut(hf.2$hist.resid,
                                  cuts), ordered=T)
  wtable2<-table(residual.magnitude2, cut(hf.1$jday,q.test))
  wtable4<-table(residual.magnitude4, cut(hf.2$jday,q.test))
  
  yrange<-range(c(as.numeric(wtable2),as.numeric(wtable4)))
  if(yrange[2]>=15){
  max.count<-ceiling(max(c(as.numeric(wtable2),as.numeric(wtable4)))/5)
  ymaxh<-max.count*5 +5
  }
  if(yrange[2]<15){
    max.count<-ceiling(max(c(as.numeric(wtable2),as.numeric(wtable4)))/3)
    ymaxh<-max.count*3 +3
  }
  if (ymaxh<=15) by.int2<-3
  if (ymaxh<30 & ymaxh>15) by.int2<-5
  if (ymaxh>=30) by.int2<-10
    
  #create plot
  plot.2<-ggplot(hf.1, aes(x=hf.1$jday,
                           fill=residual.magnitude2)
                 ,environment = .e)+
    geom_histogram(breaks=q.test) +
    scale_y_continuous(expand=c(0,0),limits= c(0,ymaxh),breaks=seq(0,ymaxh,by.int2))+
    scale_x_discrete(breaks= x.ticks,labels=label.v2)+
    coord_polar() + ggtitle("Period 1, High residuals") +
    scale_fill_discrete(name=paste("Summary statistics:\nmu=", bquote(.(hf.1.c[1])), "\nrho=", bquote(.(hf.1.c[2])),
                                   "\nkappa=", bquote(.(hf.1.c[3])),
                                   
                                   "\n\n\n\nResidual magnitude"))+
    xlab("Ordinal day")
  
  plot(plot.2)
  
  plot.3<-ggplot(lf.2, aes(x=lf.2$jday,fill=t3
  ), environment = .e)+
    geom_histogram(breaks=q.test) +
    scale_y_continuous(expand=c(0,0),limits= c(0,ymaxl),breaks=seq(0,ymaxl,by.int))+
    scale_x_discrete(breaks= x.ticks,labels=label.v2)+
    coord_polar() + ggtitle("Period 2, Low historical residuals")+
    scale_fill_discrete(name=paste("Summary statistics:\nmu=", bquote(.(lf.2.c[1])), "\nrho=", bquote(.(lf.2.c[2])),
                                   "\nkappa=", bquote(.(lf.2.c[3])),
                                   
                                   "\n\n\n\nResidual magnitude"))+
    xlab("Ordinal day")
  plot(plot.3)
  
  
  #high flow p2
  
  plot.4<-ggplot(hf.2, aes(x=hf.2$jday,fill=residual.magnitude4), environment = .e)+
    geom_histogram(breaks=q.test) +
    scale_y_continuous(expand=c(0,0),limits= c(0,ymaxh),breaks=seq(0,ymaxh,by.int2))+
    scale_x_discrete(breaks= x.ticks,labels=label.v2)+
    coord_polar() +ggtitle("Period 2, High historical residuals") +
    scale_fill_discrete(name=paste("Summary statistics:\nmu=", bquote(.(hf.2.c[1])), "\nrho=", bquote(.(hf.2.c[2])),
                                   "\nkappa=", bquote(.(hf.2.c[3])),
                                   
                                   "\n\n\n\nResidual magnitude"))+
    xlab("Ordinal day") #+
  
  
  plot(plot.4)
}




decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }

}


