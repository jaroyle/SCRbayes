\name{telemetry}
\alias{telemetry}
\docType{data}
\title{
  SCR data set that contains a sample of individuals with telemetry data
  %%   ~~ data name/kind ... ~~
}
\description{
  Simulated data from the sim.data() function
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(telemetry)}
\format{
   The format is:
 List of 7
 $ Y         : int [1:53, 1:12, 1:120] 0 0 0 0 0 0 0 0 0 1 ...

 $ MASK      : num [1:120, 1:12] 1 1 1 1 1 1 1 1 1 1 ...

 $ traplocs  : num [1:120, 1:2] 125 126 125 124 125 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:120] "1" "2" "3" "4" ...
  .. ..$ : chr [1:2] "x" "y"

 $ Xss       : num [1:947, 1] -0.566 -0.292 -0.174 0.607 -0.21 ...

 $ Ytel      : num [1:1200, 1:3] 126 126 126 126 126 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:3] "" "" "i"
  ..- attr(*, "sex")= int [1:12] 0 1 0 1 0 0 1 1 1 1 ...

 $ Xsex      : int [1:53] 0 1 1 0 1 1 1 0 0 0 ...

 $ statespace: num [1:947, 1:3] 122 122 122 123 123 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:947] "1" "2" "3" "4" ...
  .. ..$ : chr [1:3] "X_Coord" "Y_Coord" "tighab"
}
\details{
This is a list having the required objects to fit a barebones SCR model
that includes auxilliary data from a sample of telemetered
individuals. The elements of the list Y, MASK, traplocs and statespace are standard
SCR data. In addition, the following objects are provided:

Xss: a state-space covariate on "density". The data were simulated with
a strong effect of Xss determining the home range centers.

Ytel: a matrix with 3 columns having telemetry observations for a sample
of 12 individuals. The format is: x-coord, y-coord, individual_id. Ytel
has an attribute "sex" which is the sex of the 12 telemetered
individuals.

Xsex: the sex of the sample of individuals encountered by the SCR study
(not including telemetered individuals).

The data set was generated such that the telemetry and SCR samples of
individuals ARE NOT RECONCILED and assumed to be independent. 

%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Data come from sim.data()
  %%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{

Sollmann, R., Gardner, B., Chandler, R. B., Shindle, D. B., Onorato,
D. P., Royle, J. A., & O'Connell, A. F. (2013). Using multiple data
sources provides density estimates for endangered Florida
panther. Journal of Applied Ecology.

Royle, J. A., Chandler, R. B., Sun, C. C., & Fuller,
A. K. (2013). Integrating resource selection information with spatial
capture-recapture. Methods in Ecology and Evolution.

}
\examples{


### Here is how the data were simulated
data(tigerdata)
ssgrid<- tigerdata$grid900
ssgrid<- ssgrid/5000
ssgrid[,3]<- rep(1,nrow(ssgrid))
## Simulated a spatial covariate
## This is very inefficient and/or will crash your computer if ssgrid has 
##    a lot of rows.
set.seed(2014)
D<-e2dist(ssgrid,ssgrid)
V<- exp(-D/2)
x<-t(chol(V))\%*\%rnorm(nrow(ssgrid))


## Plot the covariate
par(mar=c(3,3,3,6))
spatial.plot(ssgrid,x,cx=2,add=FALSE)

traplocs<- tigerdata$tigerdata.traplocs/5000
K<- 12

## Use sim.data to simulate some data. 
telemetry<-sim.data(N = 200, sigma = c(0.35,.5), loglam0 = log(0.5), K=K, 
   statespace = ssgrid, traplocs = as.matrix(traplocs), Xss = x, alpha1 = 1, 
   coord.scale = 1,Ntel=12,nfixes=100) 

### Data can be loaded instead by doing this:
data(telemetry)

statespace<- telemetry$statespace
traplocs<-telemetry$traplocs
y<- array3d2SCR.fn(telemetry$Y)
x<- telemetry$Xss

## Reformat the "captures" data into the standard encounter data file
##    (EDF) format

edf <-cbind(session = rep(1,nrow(y)), individual=y[,2],
   occasion=y[,3],   trapid=y[,1])

## No animals were removed (dead) so "alive" is a matrix of 1's

alive=matrix(1,nrow=length(unique(edf[,"individual"])),ncol=1)

trapfile<-cbind(1:nrow(traplocs),traplocs,matrix(1,nrow=nrow(traplocs),ncol=K))

##
## Create the scrdata file:
##
scrobj<- scrData(trapfile, edf, statespace, alive = NULL, 
                 Xd = x, Ytel=telemetry$Ytel,
                 Xsex=telemetry$Xsex)

##
## Now run SCRh.fn and fit the statespace covariate on density by 
## specifying Mss=1 and indicate telemetry data is to be used by
## specifying Mtel=1. Also Msexsigma=1 says fit a model with
## sex-specific "sigma"

test<-SCRh.fn(scrobj,ni=400, burn=20, skip=1,nz=200,theta=1,
      Msigma=1, Mb=0, Msex=0, Msexsigma=1,Mss=1,Mtel=1,coord.scale=1,
      thinstatespace=1)


## maybe str(telemetry) ; plot(telemetry) ...
}
\keyword{datasets}
