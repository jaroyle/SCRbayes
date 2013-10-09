plot.scrdata<-function(obj){
if(   !any(class(obj)!="scrdata"))
{
print("not scrobj")
return(NULL)
}
coords<-as.matrix(obj$traps[,2:3])

plot(coords,pch=3)

nind<-max(obj$captures[,"individual"])
s<-matrix(NA,nrow=nind,ncol=2)
for(i in 1:nind){
trps<- obj$captures[,"trapid"][obj$captures[,"individual"]==i]
trps<- coords[trps,]
trps<-matrix(trps,ncol=2,byrow=FALSE)
s[i,]<-apply(trps,2,mean)

}
points(s,pch=20,col="red")

}