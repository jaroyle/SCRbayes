SCRdensity <-function(obj,nx=30,ny=30,Xl=NULL,Xu=NULL,Yl=NULL,Yu=NULL,
scalein=1,scaleout=1,ncolors=10, opt.ss = FALSE){
## scalein = size of pixel in area you care about
##
##
### Jan 29 2014 -- Andy dislikes the scalein and scaleout business
### trying to do something about that
##
##
##  default units is "bears per pixel" if opt.ss = TRUE



if(any(class(obj)=="scrfit")){

if (!opt.ss){
S<-obj$Sout
ss<-obj$statespace
z<-obj$zout

Sxout<-Syout<- matrix(NA,nrow=nrow(obj$Sout),ncol=ncol(obj$Sout))
for(i in 1:nrow(obj$Sout)){
tmp<- obj$Sout[i,]

Sxout[i,]<-obj$statespace[tmp,1]
Syout[i,]<-obj$statespace[tmp,2]
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
cat("sum: ",sum(Dn),fill=TRUE)
par(mar=c(3,3,3,6))
image(xg,yg,Dn,col=terrain.colors(ncolors))
image.scale(Dn,col=terrain.colors(ncolors))
box()
return(list(grid=cbind(xg,yg),Dn=Dn))

}

if (opt.ss){
# Find the centers for lions that are included in population
centers= (obj$Sout*obj$zout)
# gridchain stores the number of lions at each node of the statespace
gridchain = rep(NA,nrow(obj$statespace))

# For each point in statespace
for (j in 1:nrow(obj$statespace)){
# Sum the number of activity centers at that point
gridchain[j] = sum(centers==j)

}

n.animals=gridchain/nrow(centers) # Avg Number of Animals at each point of statespace across MCMC samples. This is density per statespace pixel

# Extract x- & y-ccordinates of statespace for image
ss.x <- sort(unique(obj$statespace[,1]))
ss.y <- sort(unique(obj$statespace[,2]))
# Find Resolution of statespace
ss.res = c(as.numeric(names(table(diff(ss.x)))),as.numeric(names(table(diff(ss.y)))))
# Calculate density scaling factor
####out.res = scaleout/prod(ss.res)*scalein
# Scale density
Dn = n.animals*scaleout  ####*out.res

cat("sum: ",sum(Dn),fill=TRUE)

ssp<-
function (x, y, add = TRUE, cx = 1)
{
    nc <- as.numeric(cut(y, 20))
    if (!add)
        plot(x, pch = " ")
    points(x, pch = 15, col = terrain.colors(20)[nc], cex = cx)
    image.scale(y, col = terrain.colors(20))
}


### Should use raster package function rasterFromXYZ here I think.

ssp(obj$statespace[,1:2][obj$statespace[,3]==1,],Dn,add=FALSE,cx=1)
###title("local density (bears / 1000)")

if(1==2){
## This code doesn't make sense if the state-space is not regular
    # Make density into a matrix
Dn = matrix(Dn, nrow = length(ss.y),ncol=length(ss.x),byrow=F)
# Plot density as before
par(mar=c(3,3,3,6))
image(ss.x,ss.y,t(Dn),col=terrain.colors(ncolors))
image.scale(Dn,col=terrain.colors(ncolors))
##image.scale(out.dens,col=terrain.colors(ncolors))
box()
}
#Return a data.frame with x-coords and y-coords and corresponding density.
#Notes that this output is slightly different from that produced above since
#the density and coordinates are of the same dimension.
#Maybe conceptually easier for user to handle? Instead of the unique x and y values.
#return(data.frame(X_coord=obj$statespace$X_coord, Y_coord=obj$statespace$Y_coord,Dn=c((Dn))))
return(cbind(obj$statespace,Dn))

}

}
else{
###
### # below lines of code for SCRbook examples

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

