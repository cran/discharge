decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
  
}


sigmaeventsplot<-function(x) {
  if(class(x)!="streamflow") stop ("Input should be streamflow object")
  flowdata<-fourierAnalysis(x)$signal
  sigma.lf<-sigmaLowFlows(flowdata, 10)
  sigma.hf<-sigmaHighFlows(flowdata, 10)
  

  .e<-environment()
  lf.1<-sigma.lf$onesigma.events
  lf.1$jday[which(lf.1$jday==366)]<-365
  hf.1<-sigma.hf$onesigma.events
  hf.1$jday[which(hf.1$jday==366)]<-365

  o<-list(sigma.hf=sigma.hf, sigma.lf=sigma.lf)
  circ<-round(circ.s(o),3)
  lf.1.c<-circ[2,2:4]
  hf.1.c<-circ[1,2:4]
  
  
  ################################
  # these define breaks for x and where to put ticks on plot
  q.test<-seq(0,365,14)
  q.test[length(q.test)]<-366
  x.ticks<-c(0,31,59,90,120,151,181,212,243,273,304,334)
  label.v2<-c("J","F","M","A","M","J","J","A","S","O","N","D")
  
  #low flow
  rres<-c(range(lf.1$resid.sig)[1], 
          range(lf.1$resid.sig)[2])
  cuts<-getbins(rres)
  
  residual.magnitude1<-factor(cut(lf.1$resid.sig,
                                  cuts),ordered=T)
  t1<-reorder(residual.magnitude1, -as.numeric(residual.magnitude1))
  
  # to set y range for low flow
  wtable1<-table(residual.magnitude1, cut(lf.1$jday,q.test))
  yrange<-range(as.numeric(wtable1))
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
  

  plot.1<-ggplot(lf.1, aes(x=lf.1$jday,fill=t1),environment = .e)+
    geom_histogram(breaks = q.test) + coord_polar() +
    scale_y_continuous(expand=c(0,0),limits= c(0,ymaxl)
                       ,breaks=seq(0,ymaxl,by.int))+
    ggtitle("Period 1, Low residuals") +
    scale_x_discrete(breaks= x.ticks,labels=label.v2)+
    scale_fill_discrete(name=paste("Summary statistics:\nmu=", bquote(.(lf.1.c[1])), 
      "\nrho=", bquote(.(lf.1.c[2])),
      "\nkappa=", bquote(.(lf.1.c[3])),
      "\n\n\n\nResidual magnitude"))+
    xlab("Ordinal day")
  par(ask=T)
  plot(plot.1)
  
  
  
  
  #low flow
  rres<-range(hf.1$resid.sig)
  cuts<-getbins(cuts)
  
  residual.magnitude2<-factor(cut(hf.1$resid.sig, 
                                  cuts),ordered=T)
  
  # to set y range for high flow
  wtable2<-table(residual.magnitude2, cut(hf.1$jday,q.test))
  
  yrange<-range(as.numeric(wtable2))
  if(yrange[2]>=15){
    max.count<-ceiling(max(c(as.numeric(wtable2),as.numeric(wtable1)))/5)
    ymaxh<-max.count*5 +5
  }
  if(yrange[2]<15){
    max.count<-ceiling(max(c(as.numeric(wtable2),as.numeric(wtable1)))/3)
    ymaxh<-max.count*3 +3
  }
  if (ymaxh<=15) by.int2<-3
  if (ymaxh<30 & ymaxh>15) by.int2<-5
  if (ymaxh>=30) by.int2<-10
  
  plot.2<-ggplot(hf.1, aes(x="jday",fill=residual.magnitude2), 
                 environment = .e)+
    geom_histogram(breaks = q.test) + coord_polar() +
    scale_y_continuous(expand=c(0,0),limits= c(0,ymaxh)
                       ,breaks=seq(0,ymaxh,by.int2))+
    ggtitle("High residuals") +
    scale_x_discrete(breaks= x.ticks,labels=label.v2)+
    scale_fill_discrete(name=paste("Summary statistics:\nmu=", bquote(.(hf.1.c[1])), 
                                   "\nrho=", bquote(.(hf.1.c[2])),
                                   "\nkappa=", bquote(.(hf.1.c[3])),
                                   "\n\n\n\nResidual magnitude"))+
    xlab("Ordinal day")
  
  plot(plot.2)

}