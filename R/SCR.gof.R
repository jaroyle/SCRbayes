SCR.gof <-
function(out,nx=20,ny=20,Xl=NULL,Xu=NULL,Yl=NULL,Yu=NULL){
## works with SCRf.fn output
S<-out$mcmchist
G<-out$G
Sxout<-Syout<-matrix(NA,nrow=nrow(S),ncol=ncol(S))
for(i in 1:nrow(S)){
Sxout[i,]<- G[,1][S[i,]]
Syout[i,]<-G[,2][S[i,]]
}
z<-out$zout
niter<-nrow(z)

if(is.null(Xl)){
Xl<-min(Sxout)*.999
Xu<-max(Sxout)*1.001
Yl<-min(Syout)*.999
Yu<-max(Syout)*1.001
}
xg<-seq(Xl,Xu,,nx)
yg<-seq(Yl,Yu,,ny)

Sxout2<-cut(Sxout[z==1],breaks=xg)
Syout2<-cut(Syout[z==1],breaks=yg)

Dn<-table(Sxout2,Syout2)/niter

image(xg,yg,Dn,col=terrain.colors(10))
image.scale(Dn,col=terrain.colors(10))

stat<-statsim<-rep(NA,niter)
for(i in 1:niter){
Dn<- table(cut(Sxout[i,][z[i,]==1],breaks=xg),cut(Syout[i,][z[i,]==1],breaks=yg))
Dnv<-Dn[1:length(Dn)]
stat[i]<-(length(Dnv)-1)*(var(Dnv)/mean(Dnv))

Sxsim<-sample(G[,1],sum(z[i,]),replace=TRUE)
Sysim<-sample(G[,2],sum(z[i,]),replace=TRUE)

Dnsim<- table(cut(Sxsim,breaks=xg),cut(Sysim,breaks=yg))
Dnsimv<-Dnsim[1:length(Dnsim)]
statsim[i]<- (length(Dnsimv)-1)*(var(Dnsimv)/mean(Dnsimv))

}


out<-cbind(data=stat,newdata=statsim)

cat("P-value: ", mean(out[,1]>out[,2]),fill=TRUE)

invisible(out)

}
