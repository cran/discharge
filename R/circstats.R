library("CircStats") 

circ.s<-function(x) {
  if(class(x)=="compflows"){
    p1.h<-cbind(circ.summary(x$p1.hevents$jday*(360/366)*pi/180), est.kappa(x$p1.hevents$jday*(360/366)*pi/180))
    colnames(p1.h)<-c("n", "mu", "rho","kappa")
    p1.l<-cbind(circ.summary(x$p1.levents$jday*(360/366)*pi/180), est.kappa(x$p1.levents$jday*(360/366)*pi/180))
    p2.h<-cbind(circ.summary(x$p2.hevents$jday*(360/366)*pi/180), est.kappa(x$p2.hevents$jday*(360/366)*pi/180))
    p2.l<-cbind(circ.summary(x$p2.levents$jday*(360/366)*pi/180), est.kappa(x$p2.levents$jday*(360/366)*pi/180))
    
    colnames(p1.l)<-c("n", "mu", "rho","kappa")
    colnames(p2.h)<-c("n", "mu", "rho","kappa")
    colnames(p2.l)<-c("n", "mu", "rho","kappa")
    
    circstats<-as.data.frame(rbind(p1.h, p1.l, p2.h, p2.l))
    row.names(circstats)<-c("p1.high", "p1.low", "p2.high", "p2.low")
  }
  
  if(class(x)!="compflows"){
    x1<-x[[1]]
    x2<-x[[2]]
    p1.h<-cbind(circ.summary(x1$onesigma.events$jday*(360/366)*pi/180), 
                est.kappa(x1$onesigma.events$jday*(360/366)*pi/180))
    colnames(p1.h)<-c("n", "mu", "rho","kappa")
    p1.l<-cbind(circ.summary(x2$onesigma.events$jday*(360/366)*pi/180), 
                est.kappa(x2$onesigma.events$jday*(360/366)*pi/180))
    
    colnames(p1.l)<-c("n", "mu", "rho","kappa")
    
    
    circstats<-as.data.frame(rbind(p1.h, p1.l))
    row.names(circstats)<-c("p1.high", "p1.low")
  }
  
  circstats
}