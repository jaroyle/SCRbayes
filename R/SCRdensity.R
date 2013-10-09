SCRdensity <-
function(obj,nx=30,ny=30,Xl=NULL,Xu=NULL,Yl=NULL,Yu=NULL,
scalein=1,scaleout=1000000*100,ncolors=10){
## 1000000 scaleout puts density in units per km^2
## multiplied by 100 is units per 100 km^2
## assumes units input are meters.

if(any(class(obj)=="scrfit")){
S<-obj$Sout
ss<-obj$statespace
z<-obj$zout

Sxout<-Syout<- matrix(NA,nrow=nrow(obj$Sout),ncol=ncol(obj$Sout))
for(i in 1:nrow(obj$Sout)){
tmp<- obj$Sout[i,]

Sxout[i,]<-obj$statespace[tmp,1]
Syout[i,]<-obj$statespace[tmp,2]
}

}
else{
# below lines of code for SCRbook examples
Sxout<-obj$Sx
Syout<-obj$Sy
z<-obj$z
}

niter<-nrow(z)

if(is.null(Xl)){
Xl<-min(Sxout)*.999
Xu<-max(Sxout)*1.001
Yl<-min(Syout)*.999
Yu<-max(Syout)*1.001
}
xg<-seq(Xl,Xu,,nx)
yg<-seq(Yl,Yu,,ny)

Sxout<-cut(Sxout[z==1],breaks=xg)
Syout<-cut(Syout[z==1],breaks=yg)
Dn<-table(Sxout,Syout)/niter  # Dn = avg # guy (posterior)
area<-  (yg[2]-yg[1])*(xg[2]-xg[1])*scalein # this is in sq km now

Dn<- (Dn/area)*scaleout   # now wolverines per 100 km
cat("mean: ",mean(Dn),fill=TRUE)
par(mar=c(3,3,3,6))
image(xg,yg,Dn,col=terrain.colors(ncolors))
image.scale(Dn,col=terrain.colors(ncolors))
box()
return(list(grid=cbind(xg,yg),Dn=Dn))
}
