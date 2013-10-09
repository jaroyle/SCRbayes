\name{SCRdensity}
\alias{SCRdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
makes a density plot from object of class "scrfit" %%  ~~function to do ... ~~
}
\description{
This makes a basic density plot given MCMC output from the main fitting
functions of SCRbayes (e.g., SCRh.fn, SCRi.fn)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRdensity(obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL, Yu = NULL, scalein = 1, scaleout = 100000*100, ncolors = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
should be of class scrfit  %%     ~~Describe \code{obj} here~~
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{Xl}{
%%     ~~Describe \code{Xl} here~~
}
  \item{Xu}{
%%     ~~Describe \code{Xu} here~~
}
  \item{Yl}{
%%     ~~Describe \code{Yl} here~~
}
  \item{Yu}{
%%     ~~Describe \code{Yu} here~~
}
  \item{scalein}{
%%     ~~Describe \code{scalein} here~~
}
  \item{scaleout}{
%%     ~~Describe \code{scaleout} here~~
}
  \item{ncolors}{
%%     ~~Describe \code{ncolors} here~~
}
}
\details{
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
  Modified from the SCR book Royle et al. (2013)
  %% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(wolverine)

## make the encounter data file (EDF)
EDF<-wolverine$wcaps
dimnames(EDF)<-list(NULL,c("session","individual","occasion","trapid"))
## grab the trap deployment file (TDF)
TDF<-wolverine$wtraps

## bundle these into an "scrdata" object along with the state-space grid
grid<-cbind(wolverine$grid2,rep(1,nrow(wolverine$grid2)))
wolv.scrdata<-scrData(TDF,EDF,grid)

## we will fit a model with sex-specific parameters, so grab the sex variable
wsex<-wolverine$wsex

## now fit a model using the current development version of
##   the fitting function
test<-SCRh.fn( wolv.scrdata,ni= 1200, burn=200, skip=2,nz=100,theta=.5,
Msigma=1, Mb=0, Msex=1, Msexsigma=1, 
coord.scale=5000, Xeff=NULL, Xsex=wsex, thinstatespace=4)

## make a density plot
SCRdensity(test, scalein=1, scaleout=1000000*100, nx=20,ny=20)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
