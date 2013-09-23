\name{make.statespace}
\alias{make.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a state-space grid for use by SCRbayes functions%%  ~~function to do ... ~~
}
\description{
This function will make a state-space grid given a set of coordinates
that define the trap locations (or points representative of potential
traps).

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.statespace(ll = NA, buffer = 0.1, minx = NA, maxx = NA, miny = NA, maxy = NA, nx = 20, ny = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ll}{
    Coordinates of traps or similar points within the study area.
    %%     ~~Describe \code{ll} here~~
}
  \item{buffer}{
    Relative size of the buffer to use in creating the state-space. A
    value of 0 produces the minimum area rectangle around the traps.
    %%     ~~Describe \code{buffer} here~~
}
\item{minx}{
  Instead of ll one could provide the minimum and maximum x and y values
  to use. 
%%     ~~Describe \code{minx} here~~
}
  \item{maxx}{
%%     ~~Describe \code{maxx} here~~
}
  \item{miny}{
%%     ~~Describe \code{miny} here~~
}
  \item{maxy}{
%%     ~~Describe \code{maxy} here~~
}
  \item{nx}{
Number of state-space points in the x-direction. Don't make this too
large. The total size of the state-space will be nx*ny points (see next
argument). A reasonable total value of nx*ny is the expected population
size N*4. So if you expect N = 100 individuals in the state-space around
the traps, nx*ny should be around 400.  
%%     ~~Describe \code{nx} here~~
}
  \item{ny=NULL}{
Number of state-space points in the y-direction. It is recommended to
leave this to NULL in which case ny is computed relative to nx.
%%     ~~Describe \code{ny} here~~
}
}
\details{
  Lots of details here......A future version of this function will have
  the user input "average home range size" and compute everything you
  need.

  The function returns the nG x 2 state-space grid which has a number of
  arguments. "area" is the area of each grid point in the units of
  "ll". "traps" is the matrix "ll". These arguments will be used by
  other functions of the SCRbayes package. 
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle,  aroyle@usgs.gov %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{






}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
