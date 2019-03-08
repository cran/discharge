
#function to create dataset with all the important info

# it will include:
# arms, nrms, snr, daily noise color, annual noise color
# sigma lf, sigma hf, log pearson iii stats:
# q2, q10, l2, l10
# zero flow?

# function to get all info for one river
allstats<-function(file.name, river.name=NULL, file.type="txt", 
                    date.col=3, discharge.col=4, skipped.rows=28) {
  cat("Starting file ", file.name, ".\n")
  
  if (file.type=="txt") {
  file<-read.table(file.name, skip=skipped.rows, sep="\t",header=F)
  }
  if (file.type=="csv") {
  file<-read.csv(file.name, skip=skipped.rows)
  }
  if (file.type!="csv" & file.type!="txt") stop ("Unrecognized file type.  
                                                Use txt or csv.")
  k<-c(date.col, discharge.col)
  file<-file[,k]
  colnames(file)<-c("date", "discharge")
  n.file<-dim(file)[1]
  max.nas<-n.file*.35
  x<-asStreamflow(file, river.name, max.na=max.nas)
  
  signal.stats<-fourierAnalysis(x)
  hflow.stats<-sigmaHighFlows(signal.stats$signal, resid.column=10)
  lflow.stats<-sigmaLowFlows(signal.stats$signal, resid.column=10)
  logpearson.stats<-lp3Events(x)
  ann.extremes<-annualExtremes(x)
  annual.stats<-annualnoise(ann.extremes$annual.max$ldis.corrupt)
  
  sigma.hf<-as.numeric(hflow.stats$sigma.hfb)
  sigma.lf<-as.numeric(lflow.stats$sigma.lfb)
  
  q2<-as.numeric(logpearson.stats[[1]])
  q10<-as.numeric(logpearson.stats[[2]])
  l2<-as.numeric(logpearson.stats[[3]])
  l10<-as.numeric(logpearson.stats[[4]])
  
  a.rms<-as.numeric(signal.stats$rms[[1]])
  n.rms<-as.numeric(signal.stats$rms[[2]])
  snr<-as.numeric(signal.stats$rms[[4]])
  theta.d<--1*as.numeric(coefficients(signal.stats$logps.regression)[2])
  name<-as.character(x$name)
  
  theta.a<--1*as.numeric(annual.stats$reg.stats[2])
  
  
  cat("File", file.name, " successful.\n")
  
  
  out<-as.data.frame(cbind(a.rms,n.rms, 
                           snr, theta.d, theta.a,sigma.hf, 
                           sigma.lf,q2, q10, l2, l10))
  
  row.names(out)<-name
  
  return(out)
}


parameters.list<-function(x, names=NULL, file.type="txt", date.col=3, 
                          dis.col=4,skipped.rows=28) {
  
  n.files<-length(x)
  files<-x
  if( !is.null(names)) name.vec<-names
  if (is.null(names)) name.vec<-as.character(1:n.files)
  
  all.out<-NULL
  filetype<-file.type
  datecol<-date.col
  discol<-dis.col
  skipped<-skipped.rows
  for (i in 1:n.files) {
  output<-tryCatch(allstats(files[i], name.vec[i], file.type=filetype, date.col=datecol, 
                             discharge.col=discol,skipped.rows=skipped), 
                   
    error=function(ex){
    cat("Error in file ", files[i])
    print(ex)
    return(rep(NA,11))})
  all.out<-rbind(all.out, output)
  }
  all.out<-as.data.frame(all.out)
  return(all.out)

}


