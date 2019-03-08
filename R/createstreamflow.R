 # function to get series into correct date format
# and normalize discharge
# would be nice to support different formats of date input
# error if start date specified and too many NAs
# if start date not specified, could it find a good period to use?

asStreamflow<-function(x, river.name=NULL, start.date=NULL, end.date=NULL,max.na=10)
{

colnames(x)<-c("date", "discharge")
n<-dim(x)[1]

#river<-x
# find nas and na blocks
nas<-which(is.na(x[,2]))
rle.nas<-rle(is.na(x[,2]))
na.lengths<-rle.nas$lengths[which(rle.nas$values==T)]
long.na<-which(na.lengths>20)

if (length(long.na)>0){
nas.sum<-cumsum(rle.nas$length)
na.starts<-nas.sum[which(rle.nas$values==T)-1]+1
na.block.starts<-as.character(x$date[na.starts[long.na]])
#na.block.starts<-as.Date(x$date, "&Y-&m-&d")
na.info<-as.data.frame(cbind(na.starts[long.na],na.block.starts, na.lengths[long.na]))
names(na.info)<-c("index","date", "length")
}
if (length(long.na)==0) na.info<-NULL




## to keep non-numeric data from being read in as a factor ##
# these lines will produce a warning if there are NAs in the data #
x$discharge<-as.character(x$discharge)
x$discharge<-as.numeric(x$discharge)

jday<-strptime(x$date, "%Y-%m-%d")$yday+1

#get numeric dates, months, years ##
x$date<-as.Date(x$date, "%Y-%m-%d")

date.num<-as.numeric(x$date)
year<-as.numeric(format(x$date, format = "%Y"))
month<-as.numeric(format(x$date, format = "%m"))
day<-as.numeric(format(x$date, format = "%d"))

# function to find 7 day moving average
ma <- function(x,n=7){filter(x,rep(1/n,n), sides=2)}
#add normalized log discharge ###
log.discharge<-log10(x$discharge+1)
mean.discharge<-mean(log.discharge, na.rm=TRUE)
ldis.corrupt<-(log10(x$discharge+1)) / mean.discharge
smoothed.7d<-as.numeric(ma(ldis.corrupt))
x<-cbind(x, date.num, year, month, day, jday, ldis.corrupt, smoothed.7d) 



#subset data using subsetData function

if(is.null(start.date)) start.date<-x$date[1]
if(is.null(end.date)) end.date<-x$date[n]


subset.data<-subsetData(start.date, end.date, x, max.nas=max.na)
sub.data<-subset.data$sub.data #to get subsetted dataframe
n.nas<-subset.data$n.nas

stream.obj<-list(data=sub.data, n=dim(sub.data)[1], n.nas=n.nas,
                 start=as.character(start.date), 
end=as.character(end.date), name=river.name, na.info=na.info)
class(stream.obj)<-"streamflow"
return(stream.obj)

}

summary.streamflow<-function(object, ...)
{
cat(object$name, "\n")
cat(object$n, " observations", "\n")
cat("Start date:", object$start, "   End date:", object$end, "\n")
cat(object$n.nas, " missing values", "\n")
cat("NA blocks \n")
object$na.info
}

print.streamflow<-function(x, ...)
{
cat(summary(x), "\n")
cat("To view all data, use $data \n \n")
print(x$data[1:5,], ...)
}