\name{captures2alive}
\alias{captures2alive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make an alive matrix from the captures matrix
%%  ~~function to do ... ~~
}
\description{
	A function that takes the captures matrix given to scrData and the status ("Live"/"Dead") of each of those captures and makes a matrix of 0's and 1's defining whether an individual was alive or dead during each sampling period.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
captures2alive(scrCaptures, status)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrCaptures}{
  	The capture matrix given to scrData. Matrix having encounter information about each individual. Some flexibility 
is allowed in data format (more flexibility in the future).
Minimally, the matrix should have columns labeled 
"individual" "occasion" and "trapid".

"individual" = integer from 1:nind

"occasion" = integer from 1:(noccasions)

"trapid" = integer from 1:ntraps mapping to rows of "traps"

Optionally may have a column "session" which is the year or primary
period of the sample. Currently no multi-session functionality is available
but the main functions of SCRbayes accepts this data structure. 
%%     ~~Describe \code{scrCaptures} here~~
}
  \item{status}{
  	A character vector of "Live" or "Dead" describing whether each capture was a live capture or a dead recovery.
%%     ~~Describe \code{status} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A matrix describing whether each individual was available for 
sampling (1) or not (0).  Models in SCRbayes assume you know 
whether individual were removed or not, without error. 
The matrix will be "nind x noccasions". 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Joshua Goldberg, joshua.goldberg@umontana.edu
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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

# When did captures occur
SO<-c(1,2,3)
for(i in 1:nrow(captures)){
  temp<-sample(SO, captures$number[i])
   for(j in 1:length(temp)){
    captures[i,temp[j]+3]<-1
   }
 }


capture.sims<-data.frame(trapid=rep(captures[,1],3), individual=rep(captures[,2],3), occasion=c(captures[,4], captures[,5]*2, captures[,6]*3)) 
capture.sims<-subset(capture.sims, capture.sims$occasion!=0)  ###check captures, eliminate SO=0
capture.sims = capture.sims[order(capture.sims$individual),] # Order by individual
capture.sims$individual = as.numeric(as.factor(capture.sims$individual)) # Remove missing individual ids

# Simulate capture status for each individual
status = sample(c("Live","Dead"),length(unique(capture.sims$individual)),replace=T)
# Assign dead captures to final sample for each individual
status.sims = rep(NA, nrow(capture.sims))
for (i in 1:length(status)){
	if (status[i]=="Live"){
		status.sims[which(capture.sims$individual==unique(capture.sims$individual)[i])]="Live"
	}
	if (status[i]=="Dead"){
		status.sims[which(capture.sims$individual==unique(capture.sims$individual)[i])]=c(rep("Live",sum(capture.sims$individual==unique(capture.sims$individual)[i])-1),"Dead")
	}
}

# Make alive matrix given captures and status
Alive.sims = captures2alive(capture.sims,status.sims)




Xsex<-rbinom(length(unique(capture.sims$individual)),1,0.5)


scrobj_harv<-scrData(traps=traps.sims, captures=capture.sims, statespace=ssanimal.sims, alive=Alive.sims, Xd=NULL) 

sim.test1<-SCRh.fn(scrobj_harv,ni=200, burn=20,nz=1000, theta=1,
thinstatespace=2, Msex=0)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
