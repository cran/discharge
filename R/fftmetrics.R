
fftmetrics=function(x,year,candmin,candmax){

sf=fourierAnalysis(x)
x.m=sf$signal
yr=year
ann.ind=c(which(x.m$year==(yr-1) & x.m$month>8),which(x.m$year==yr),which(x.m$year==(yr+1) & x.m$month<5))
ann.ind2=which(x.m$year==yr)

# search for extrema of signals using candmin/candmax
searchmin=which(!is.na(match(x.m$jday[ann.ind],candmin)))
searchmax=which(!is.na(match(x.m$jday[ann.ind],candmax)))
nblocks=rle(c(-10,searchmin)-c(searchmin,-10))
nblocks2=rle(c(-10,searchmax)-c(searchmax,-10))
nnblocks=length(which(nblocks$lengths>1))
nnblocks2=length(which(nblocks2$lengths>1)) # how many local mins
start.min=cumsum(nblocks$lengths)[which(nblocks$lengths>1)-1]
start.max=cumsum(nblocks2$lengths)[which(nblocks2$lengths>1)-1]

m=rep(NA,nnblocks)
mx=rep(NA,nnblocks)
start.min=c(start.min,length(searchmin))
start.max=c(start.max,length(searchmax))

for(i in 1:nnblocks){
  temp=searchmin[start.min[i]:(start.min[i+1]-1)]
  m[i]=temp[which.min(x.m[ann.ind,11][temp])]}
for(i in 1:nnblocks2){
  temp=searchmax[start.max[i]:(start.max[i+1]-1)]
  mx[i]=temp[which.max(x.m[ann.ind,11][temp])]}

if(length(mx)>1){
  lowflow.window=ann.ind[mx[1:2]]}
if(length(m)>1){
  highflow.window=ann.ind[m[1:2]]}
hf.wind=highflow.window[1]:highflow.window[2]
lf.wind=lowflow.window[1]:lowflow.window[2]
ref.point=lf.wind[which(x.m[lowflow.window,4]==yr)][1] # first local max in the specified year
plot.ref.point=mx[which(x.m[ann.ind[mx],4]==yr)]

# find min and max events
min.sam=which.min(x.m[ann.ind2,10])
max.sam=which.max(x.m[ann.ind2,10])
# return columns corresponding to date, ldis, resid, predicted, jday
index=c(ann.ind2[min.sam],ann.ind2[max.sam])
sam=as.data.frame(cbind(x.m[index,1],x.m[index,7], x.m[index,c(10,11)],index,
                        index-ref.point))
names(sam)[c(1,2,6)]=c("date","jday","time ref.pt")
rownames(sam)=c("minSAM:", "maxSAM:")


# to get duration of longest sequence
runs=rle(as.numeric(x.m[ann.ind,10]>0))
longest.high=runs$length[runs$values==1][which.max(runs$length[runs$values==1])]
longest.low=runs$length[runs$values==0][which.max(runs$length[runs$values==0])]
start.highrun=cumsum(runs$length)[runs$values==0][which.max(runs$length[runs$values==1])-as.numeric(runs$values[1]==1)]+1
start.lowrun=cumsum(runs$length)[runs$values==1][which.max(runs$length[runs$values==0])-as.numeric(runs$values[1]==0)]+1

q=rbind(x.m[ann.ind[start.lowrun],c(1,7,10,11)],
x.m[ann.ind[start.highrun],c(1,7,10,11)])
rownames(q)=c("lowflow event","highflow event")
q=cbind(q,c(longest.low,longest.high),c(ann.ind[start.lowrun],ann.ind[start.highrun]),-ref.point+c(ann.ind[start.lowrun],ann.ind[start.highrun]))
names(q)[c(1,5,6,7)]=c("startdate","runlength","index","time.ref.pt")


# net and relative auc
net.auc=sum(as.numeric(lapply(x.m[hf.wind,10],max,0)))+
  sum(as.numeric(lapply(x.m[lf.wind,10],min,0)))
rel.auc.low=log10(-sum(as.numeric(lapply(x.m[lf.wind,10],max,0)))/sum(as.numeric(lapply(x.m[lf.wind,10],min,0))))
rel.auc.high=log10(-sum(as.numeric(lapply(x.m[hf.wind,10],max,0)))/sum(as.numeric(lapply(x.m[hf.wind,10],min,0))))

auc=data.frame(net.auc=net.auc,rel.auc.low=rel.auc.low,
           rel.auc.high=rel.auc.high)

# last step: prep plot
days=x.m[ann.ind,7]
d.seq=seq(1,595,45)+1
plot(x.m[ann.ind,8],cex=.75,main=paste(x$name, ": ",yr),
     ylab="log normalized discharge",xlab="ordinal day",xaxt="n")
axis(1,at=d.seq,as.character(days[d.seq]))
points(x.m[ann.ind,11],type="l",lwd=3,col="red")
abline(v=which(x.m[ann.ind,4]==yr)[1],lty=2,col="grey")
abline(v=which(x.m[ann.ind,4]==yr+1)[1]-1,lty=2,col="grey")
#abline(v=plot.ref.point,col="red",lty=4)
points(plot.ref.point,x.m[ref.point,11],col="purple",cex=1.25,pch=8)
points(which(x.m[ann.ind,4]==yr)[min.sam],x.m[ann.ind2[min.sam],8],col="orange",cex=1.25,pch=18)
points(which(x.m[ann.ind,4]==yr)[max.sam],x.m[ann.ind2[max.sam],8],col="blue",cex=1.25,pch=18)
legend("topright",col=c("orange",'blue'),pch=c(18,18),
       legend=c(paste(sam[1,1]),paste(sam[2,1])))
text(which(x.m[ann.ind,4]==yr)[1]+1,par("yaxp")[2]*.85,pos=1,labels=paste(yr),cex=.75)
text(which(x.m[ann.ind,4]==yr+1)[1]+1,par("yaxp")[2]*.85,pos=1,labels=paste(yr+1),cex=.75)

l=list(sam=sam,refpoint=x.m[ref.point,c(1,2,7,10,11)],events=q,auc=auc,
       noise.color=-as.numeric(coefficients(sf$logps.regression)[2]))
return(l)
}