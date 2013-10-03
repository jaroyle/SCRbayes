\name{SCRh.fn}
\alias{SCRh.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
conducts spatial capture recapture analysis 
}
\description{
Estimates density of animals from mark-recapture surveys where individual animals were unique identfied
 and spatial locations of captured animals were recorded.Estimation is conducted in a Bayesian framework.  
}
\usage{
SCRh.fn(scrobj,ni = 1100, burn = 100, skip = 2, nz = 200,theta=2, Msigma= 1, Mb = 0, Msex = 0, Msexsigma = 0,
Xeff = NULL, Xsex = NULL, coord.scale = 5000, thinstatespace = 1, maxNN = 20, dumprate = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{scrobj}{

A list of type "scrobj" which has 5 elements "traps",
  "captures", "statespace", "alive", "Xd"

traps: traps is a matrix containing a seq of numbers (1,2,3....j) where j= the
number of traps in column 1.  Traps can be cameras, hair snares, or
centerpoints of grid locations.  In columns two and three of the trap
matrix are the X and Y or longitude and latitude of the trap locations,
and the survey dates in columns 4 to column 3+k where k=the number of
survey dates. Each row of the matrix will contain either a 1 if the trap
was operational on the survey date or a 0 if the trap was not operational.


XXXX captures -- NO, andy changed this on 7/29 XXXXXXX
captures: a matrix
where the number of rows=the number of captures.  Each row contains 3
columns labeled LOC_ID, ANIMAL_ID, and SO; where LOC_ID =the number of
the trap where the animal was captured, ANIMAL_ID = the animal's unique
number (must be a sequence 1,2,,,,number of animals), and SO=the
sampling occasion on which the capture took place 

statespace: contains the coordinate of potential centers of activity
including those from outside the trapping area, the statespace is a
matrix where the X or longitudinal coordinate is in the first column and
the Y or latitudinal coordinate is in the second column, and HABITAT= to
1 or 0 is the third column.  The third column indicates whether the
point defined by the X,Y locations is considered habitat (normally)
habitat=1; (see function make.statespace) 

alive: a matrix of (nind x K) where alive[i,k] = 1 if the individual is
thought to be alive and alive[i,k]=0 if it is known dead.

Xd = a single state-space covariate. A vector of length =
nrows(statespace).

}

\item{ni}{ ni=number of iterations }

\item{burn}{ burn=the number of iterations that are considered "burn-in"
  and will be discarded from the final estimation }

\item{skip}{
  skip=the
number of iterations in a series that will be retained.  Skip=3 will
retain every third iteration of the MCMC chain.
}
\item{nz}{ 
nz=the
number of all zero encounter histories; the number of animals to augment
the data with; experimentation with this number is recommended - the
distribution of the parameter estimate of animal density should not be
skewed or "truncated" by the number of augmented animals
 }

 \item{theta}{
   the exponent parameter of the distance function. theta = 1
  is exponential; theta=2 is Gaussian; theta = NA then theta is
  estimated.
  }

\item{Msigma}{
Msigma is an indicator variable that =1 if you are estimating the
spatial model, normally Msigma=1
}
  \item{Mb}{
Mb is an indicator variable that=1 if behavioral effects on capture histories are expected.  For example trap-happiness or trap shyness.  Setting Mb=1 is analagous to running
model Mb described by Otis et al. 1978
}
  \item{Msex}{
Msex is an indicator variable that=1 if detection probability is to be estimated seperately by sex, if Msex=1, you must provide a vector Xsex (see Xsex for more details)
}
  \item{Msexsigma}{
Msex sigma is an indicatore variable that = 1 if sigma is to be estimated seperately by sex, if Msexsigma =1 you must provide a vector Xsex
}

  \item{Xeff}{ Xeff=j by k matrix where k=the number of surveys and
j=the number of traps.  For each row indicate the amount of survey
effort conducted in that trap location during survey k, if effort was
equal across all surveyed areas Xeff=NULL; If a trap x occasion was not
sampled then any value can be used for effort because those elements of
Xeff will not be extracted for the likelihood calculations. (the file
"traps" contains information about which traps x occasions were sampled,
see above).  }

\item{Xsex}{ Xsex= a vector of length n where n= the
number of animals;Xsex[n] is either a zero or 1 indicating the sex of
animal [n] in the capture file, must be provided if Msex or Msexsigma=1
} 

\item{coord.scale}{ coord.scale scales the coordinate system by the
input value.  The default is 5000 and will scale all X and Y locations
by 5000 }

\item{area.per.pixel=1}{ Area of a state-space grid point. Default is 1
  unit of area. This is used in the conversion of population size N to
  density D. Using this value, the total area of the state-space is
  computed and D = N/totalarea is a derived parameter. Posterior
  summaries of D will be reported by print() and other summary
  functions. 
  If thinstatespace is specified then the area.per.pixel
  is adjusted.

  If coordinates are UTM or any other scale but you wish to have density
  reported in some other units, then provide  area.per.pixel in the
  relevant units. e.g., if a grid cell is 1000 x 1000 meters and you
  wish density to be in units of "per km^2" then specify area.per.pixel
  = 1. If you want density in units "per 100 km^2" then specify
  area.per.pixel = 0.01. 
}

\item{thinstatespace}{ thinstate allows the user to retain
every yth point in the state space, fewer points can add with
computational time, but two few point will bias density estimates }

\item{maxNN}{ used to determine the "neighborhood" for each grid cell;
getNN() function in SCRh.fn; information is used in updating the
activity centers in the MCMC algorithm }

\item{dumprate}{ number of
  iterations to run before writing a file of results }

\item{call}{ executed function call }
}

\details{ User is
advised to read referenced literature, and experimentation with multiple
values to determine the appropriate nz, and the size and density of the
statespace.  SCRh.fn will provide a summary at the beginning telling the
user how many animals were captured and how many were recaptured.
User's should be aware of limitations of their data and determine for
themselves whether SCR methods are appropriate for their study.



Model implemented by SCRh.fn:

The basic encounter model that exists in SCRh.fn is:

 Pr(y=1) = 1-exp(-lam0*k(x,s))

where

 k(x,s) = exp(-(1/(2*sigma^2))*dist(x,s)^(2/theta))
 
Therefore, for theta = 1/2 this is an exponeital kernal and for theta = 1 this is a Gaussian kernel.  
Note that the encounter model is usually called the hazard model or complementary log-log model because:

cloglog(Pr(y=1)) = log(lam0) + log(k(x,s))
 
which is why you see things in the code like this:

lp.sigmac<- Msigma*bsigmac*(c1+c2)^theta
lpc<-  loglam0 + Mb*beta.behave*prevcap - lp.sigmac + beta1*Xeff  + Msex*beta.sex*Xsex[indid]

(lp being short for "linear predictor").

We see that 

lpc  = bunchofjunk  - func.of.sigma*(distance.squared)^theta
  

so that:

exp(lpc) = constant*exp(-func.of.sigma*(distance^2/theta))

 

 

 


}


\value{ returns an object with several features.  the function names can
be used to identify the features associated with the returned object.
features that can be accessed using $ include, out- a matrix containing
the values of the parameter estimates for all retained iterations.

G - a matrix containing the scaled coordinates of the statespace,

traplocs - a matrix containing the scaled coordinates of the traps, 


Sout - an i by n matrix containing the number corresponding to the row in the statespace that contains the coordinates of the estimated activity center for animal n in iteration i

zout - an i by n matrix indicating where the animal was estimated to be a member of the population zout=1 or not zout=0 for each i iteration of the model; zout combined with Sout can be
used to determine the locations of the animals that were included in the population.see functions spatial.plot and image.scale

}
\references{
Otis, D. L., K. P. Burnham, G. C. White, and D. R. Anderson. 1978. Statistical inference from capture
data on closed animal populations. Wildlife Monographs 62.

Royle, J. A., A. J. Magoun, B. Gardner, P. Valkenburg, and R. E. Lowell. 2011. 
Density estimation in a wolverine population using spatial capture-recapture models.  
Journal of Wildlife Management 75:604-611.

Gardner, B., J. A. Royle, and M. T. Wegan. 2009. Hierarchical models for estimating density from DNA mark–recapture studies.
 Ecology 90:1106-1115.

Thompson, C. M., J. A. Royle, J. D. Gardner. in press. A framework for inference about carnivore density from unstructured spatial sampling of scat using detector dogs.
  Journal of Wildlife Management.
}
\author{
Andy Royle and Robin Russell 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


## Load the mountain lion data
data("lions")

## Reformat the "captures" data into the standard encounter data file (EDF) format
newcaptures.lions<-cbind(
   session = rep(1,nrow(captures.lions)),
   individual=captures.lions[,2],
   occasion=captures.lions[,3],
   trapid=captures.lions[,1])

## No animals were removed (dead) so "alive" is a matrix of 1's
alive=matrix(1,nrow=length(unique(newcaptures.lions[,"individual"])),ncol=1)
## State-space covariate: just constant in this case
Xd<- rep(1,nrow(statespace.lions))

### old format
### scrobj<-list(traps=traps.lions,captures=captures.lions,statespace=statespace.lions)
scrobj<-list(
    traps=traps.lions,
    captures=newcaptures.lions,
    statespace=statespace.lions,
    alive=alive,
    Xd=Xd)

test<-SCRh.fn(scrobj,ni=220, burn=20, skip=2,nz=200,
      Msigma=1, Mb=0, Msex=0, Msexsigma=0,
      Xeff=NULL, Xsex=NULL, thinstatespace=4)

####
###
### here are 2 comprehensive examples based on simulated data. This can be
###    used to test the function and make sure its working. 
####
####

#### make a 10 by 10 grid of traps  km apart

min_y=0
max_y=100
min_x=0
max_x=100
grid.density=10

y<-seq(min_y, max_y, by=grid.density)
x<-seq(min_x, max_x, by=grid.density)

trap_locs<-expand.grid(x[1:length(x)], y[1:length(y)])
coord.scale=10
ntraps<- nrow(trap_locs)
 
traps.sims<-data.frame(TRAP_ID=seq(1,nrow(trap_locs)),x=trap_locs[,1], y=trap_locs[,2], 
    SO1=rep(1, nrow(trap_locs)), So2=rep(1, nrow(trap_locs)), So2=rep(1, nrow(trap_locs)))

buffer=10
make.ss<-function(traps, buffer, grid.density){
x<-seq(min(traps.sims$x)-buffer, max(traps.sims$x)+buffer, by=grid.density)
y<-seq(min(traps.sims$y)-buffer, max(traps.sims$y)+buffer, by=grid.density)
ss.y<-rep(y, length(x))
ss.x<-rep(x[1], length(y))
for (i in 2:length(x)) {
  ss.x<-c(ss.x,rep(x[i],length(y)))
}
ss.animal<-data.frame(X_coord=ss.x, Y_coord=ss.y, V3=1)
list(ss.animal=ss.animal)
}

##make a state space, buffer= buffer + 1/2 grid
out<-make.ss(traps.sims, buffer,1)
ssanimal.sims<-out$ss.animal

###simulate activity centers 
N<-1000
z_centers<-data.frame(x=rep(0, N), y=(rep(0,N)))
Potential_x<-seq(-10,110, by=1)
Potential_y<-seq(-10,110, by=1)

for(i in 1:N){
   z_centers[i,1]<-sample(Potential_x, 1)
   z_centers[i,2]<-sample(Potential_y,1)
}

plot(ssanimal.sims[,1], ssanimal.sims[,2])
points(traps.sims$x, traps.sims$y, col="red", bg="red", pch=21)
points(z_centers[,1],z_centers[,2], col="blue", bg="blue", pch=21)
 
######calculate probability of capture for individual i in trap j
###basics for sigma= 4,5,6,7,8,
sigma=10
d<-e2dist(z_centers,traps.sims[,2:3])
alpha0 <- -2.5
alpha1 <- 1/(2*sigma*sigma)
prob_cap<- plogis(-2.5)*exp( - alpha1*d*d)

y<-matrix(NA,nrow=N,ncol=ntraps)

for(i in 1:nrow(y)){
  y[i,]<-rbinom(ntraps,3,prob_cap[i,])
}

captures<-data.frame(trapid=col(y)[y>0],individual=row(y)[y>0],number=y[y>0])
captures$SO=rep(0, nrow(captures))
captures$SO1=rep(0, nrow(captures))
captures$SO2=rep(0, nrow(captures))

SO<-c(1,2,3)
for(i in 1:nrow(captures)){
  temp<-sample(SO, captures$number[i])
   for(j in 1:length(temp)){
    captures[i,temp[j]+3]<-1
   }
 }

#####pick h for harvest
H<-200
harvest<-sample(1:N,H)
###randomly assign when they were harvested
harvest.time<-sample(c(1,2), replace=TRUE,length(harvest))

########set captures to zero that occured after harvest 
for(i in 1:nrow(captures)){
 for(j in 1:length(harvest)){
    captures[i,5]<-ifelse(captures$individual[i]==harvest[j] & harvest.time[j]==1,0, captures[i,5])
    captures[i,6]<-ifelse(captures$individual[i]==harvest[j] & harvest.time[j]<3,0, captures[i,6])
  }
 }

###########create Alive.sims matrix
Alive.sims<-data.frame(individual=unique(captures$individual))
Alive.sims[,2:4]<-1
for(i in 1:nrow(Alive.sims)){
  for(j in 1:length(harvest)){
    Alive.sims[i,3]<-ifelse(Alive.sims$individual[i]==harvest[j] & harvest.time[j]==1,0, Alive.sims[i,3])
    Alive.sims[i,4]<-ifelse(Alive.sims$individual[i]==harvest[j] & harvest.time[j]<3,0, Alive.sims[i,4])
  }
}

capture.sims<-data.frame(trapid=rep(captures[,1],3), individual=rep(captures[,2],3), occasion=c(captures[,4], captures[,5]*2, captures[,6]*3)) 
capture.sims<-subset(capture.sims, capture.sims$occasion!=0)  ###check captures, eliminate SO=0


#######don't forget to put the alive.sims matrix in order!!!
Alive.sims<- Alive.sims[ order(Alive.sims$individual) ,]
Xsex<-rbinom(length(unique(capture.sims$individual)),1,0.5)


scrobj_harv<-scrData(traps=traps.sims, captures=capture.sims, statespace=ssanimal.sims, alive=Alive.sims, Xd=NULL) 
sim.test1<-SCRh.fn(ni=200, burn=20,scrobj_harv,nz=1000, theta=1,thinstatespace=2, Msex=0, Xsex=NULL, Xeff=NULL)



####
#### old code
####
####

N<-200
 
traps.sims<-data.frame(TRAP_ID=seq(1,nrow(trap_locs)),x=trap_locs[,1], y=trap_locs[,2], SO1=rep(1, nrow(trap_locs)), So2=rep(1, nrow(trap_locs)), So2=rep(1, nrow(trap_locs)))
ntraps<-nrow(traps.sims)

##make a state space, buffer= buffer + 1/2 grid
out<-make.ss(traps, buffer,1)

ss.animal.sims<-out$ss.animal

###simulate activity centers 
z_centers<-data.frame(x=rep(0, N), y=(rep(0,N)))
Potential_x<-seq(-10,110, by=1)
Potential_y<-seq(-10,110, by=1)

### a fake covariate here
f<- .2*(ss.animal.sims[,1] -50)/20  + .2*(ss.animal.sims[,2] - 50)/20 + -.25*((ss.animal.sims[,1] - 50)/20)*((ss.animal.sims[,2] - 50)/20)
spatial.plot(ss.animal.sims[,1:2],f,add=FALSE)

beta.den<-1
probs<-exp(beta.den*f)/sum(exp(beta.den*f))

Zmat<- ss.animal.sims[sample( 1:nrow(ss.animal.sims), N, replace=TRUE, prob=probs),1:2]

z_centers[,1] <- Zmat[,1]
z_centers[,2] <- Zmat[,2]

plot(ss.animal.sims[,1], ss.animal.sims[,2])
points(traps$x, traps$y, col="red", bg="red", pch=21)
points(z_centers[,1],z_centers[,2], col="blue", bg="blue", pch=21)

 
######calculate probability of capture for individual i in trap j

sigma=10

d<-e2dist(z_centers,traps.sims[,2:3])

alpha0 <- -2.5
alpha1 <- 1/(2*sigma*sigma)
prob_cap<- plogis(-2.5)*exp( - alpha1*d*d)

K<-3
y<-matrix(NA,nrow=N,ncol=ntraps)

for(i in 1:nrow(y)){
   y[i,]<-rbinom(ntraps,K,prob_cap[i,])
}

captures<-data.frame(trapid=col(y)[y>0],individual=row(y)[y>0],number=y[y>0])
captures$SO=rep(0, nrow(captures))
captures$SO1=rep(0, nrow(captures))
captures$SO2=rep(0, nrow(captures))

SO<-c(1,2,3)
for(i in 1:nrow(captures)){
  temp<-sample(SO, captures$number[i])
  for(j in 1:length(temp)){
     captures[i,temp[j]+3]<-1
   }
}


### No harvest in this example
Alive.sims<-matrix(1,nrow=length(unique(captures[,2])),ncol=K)

session<- rep(1,nrow(captures))

capture.sims<-data.frame(
session=session,
individual=rep(captures[,2],3), 
SO=c(captures[,4], captures[,5]*2, captures[,6]*3),
trapid=rep(captures[,1],3)
) 

capture.sims<-subset(capture.sims, capture.sims$SO!=0)  ###check captures, eliminate SO=0

### This example does NOT use "scrData"
scrobj2<-list(traps=traps.sims, captures=capture.sims, statespace=ss.animal.sims, alive=Alive.sims,Xd=f)
Xsex<-rbinom(length(unique(capture.sims$individual)),1,0.5)

time1<- unix.time(
  test2<-SCRh.fn(scrobj2,ni= 1200, burn=200,
     skip=2,nz=240,theta=.5,Msigma=1, Mb=0, Msex=0, Msexsigma=0, coord.scale=1,
     Xeff=NULL, Xsex=NULL, thinstatespace=4)
)





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line




