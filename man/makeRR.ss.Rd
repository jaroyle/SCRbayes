\\name{make.ss}
\alias{make.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
this function is a support function for SCRf.fn and creates a statespace with a user specified buffer and grid density
}
\description{
function creates a square statespace around a survey area containing x,y locations of potential activity centers for use in the function
SCRf.fn.  
}
\usage{
make.ss(traps, buffer, grid.density)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traps}{
traps=matrix containing the X,Y locations of traps laid out in a grid
}
  \item{buffer}{
the distance of the buffer to put around the traps.  Remember if your trap location is the centerpoint of a grid square, add half
the size of the grid square to your desired buffer.  (i.e. for a buffer of size 10km and grid of 5X5km squares, your buffer size 
should be 12.5km
}
  \item{ss.density}{
the density of your desired statespace on the same scale as the trapping grid. An ss.density of 2 in statespace for a grid of 5x5km
will return a statespace with a point every  2km 
}


}
\details{
function is designed for rectangular trapping arrays
}
\value{
returns a matrix with three column, column 1 and 2 are the x and y coordinates of the state space, this coordinates are the locations of the potential activity
centers, column three is by default =1 and indicates whether the coordinate is in a location considered habitat or not. 
}
\references{
Royle, J. A., K. U. Karanth, A. M. Gopalaswamy, and N. S. Kumar. 2009. Bayesian inference in camera trapping studies for a class of spatial capture–recapture models. Ecology 90:3233–3244.

Thompson, C. M., J. A. Royle, J. D. Gardner. in press. A framework for inference about carnivore density from unstructured spatial sampling of scat using detector dogs.  Journal of Wildlife Management.
}
\author{
Robin Russell
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
SCRf.fn, image.scale, e2dist
}
\examples{
###create a trapping grid
traps<-data.frame(x=rep(seq(4000,5000, by=100),6), y=sort(rep(seq(4000,5000, by=200),11)))
##make a rectangular statespace with a center of activity every 100m, a buffer of 100m, and a habitat value of 1 for each center point
ss<-make.ss(traps, 100, 100)





make.ss<-function(traps, buffer, grid.density)
{
x<-seq(min(traps$x)-buffer, max(traps$x)+buffer, by=grid.density)
y<-seq(min(traps$y)-buffer, max(traps$y)+buffer, by=grid.density)
ss.y<-rep(y, length(x))
ss.x<-rep(x[1], length(y))

for (i in 2:length(x))
{
ss.x<-c(ss.x,rep(x[i],length(y)))
}

ss.animal<-data.frame(X_coord=ss.x, Y_coord=ss.y, V3=1)
list(ss.animal=ss.animal)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
