\name{e2dist}
\alias{e2dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cc%% a support function for SCRf.fn, calculates the distances between all points in a grid and a set of locations from within the grid
}
\description{
%%  function estimates the euclidian distance between the potential activity centers in the statespace and the observed locations of the captured animals
}
\usage{
e2dist(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     x is a grid or "statespace" of uniform density
}
  \item{y}{
%%    y is a set of capture locations  from within the state space
}
}
\details{
%%  e2dist is used to determine the distance between potential activity centers in the statespace and the observed location of animals from spatial mark-recapture studies
}
\value{
%%   a matrix nrows=the number of rows in grid x and ncols=the number of capture locations.  matrix contains the euclidian distance between the point in the grid and the capture location
}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  J. Andrew Royle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

x<-seq(1,100, by=10)  ###create a grid called "x"
y<-seq(1,200, by=10)	
x.y<-rep(y, length(x))
x.x<-rep(x[1], length(y))

for (i in 2:length(x))
{
x.x<-c(x.x,rep(x[i],length(y)))
}
x<-cbind(x.x,x.y)  ### x is a 
y<-cbind(sample(x[,1],50), sample(x[,1],50))  ##y is a random set of points drawn from x

## The function is currently defined as
function (x, y) 
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
